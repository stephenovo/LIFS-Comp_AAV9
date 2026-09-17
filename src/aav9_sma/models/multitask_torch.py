"""Optional PyTorch multi-task model with per-endpoint missing-label masks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _torch_modules():
    try:
        import torch
        from torch import nn
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as error:
        raise ImportError(
            "PyTorch is optional; install the project with the 'models' extra"
        ) from error
    return torch, nn, DataLoader, TensorDataset


@dataclass(frozen=True)
class MaskedMLPResult:
    """Predictions plus reproducibility metadata from one fitted model."""

    predictions: np.ndarray
    best_epoch: int
    validation_loss: float
    training_rows: int
    validation_rows: int
    task_observations: tuple[int, ...]
    device: str


def _fit_scaler(targets: np.ndarray, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    selected = np.where(rows[:, None] & np.isfinite(targets), targets, np.nan)
    means = np.nanmean(selected, axis=0)
    scales = np.nanstd(selected, axis=0)
    scales = np.where(scales > 0, scales, 1.0)
    if not np.isfinite(means).all():
        raise ValueError("Every task needs at least one observed training label")
    return means, scales


def _train_epochs(
    features: np.ndarray,
    targets: np.ndarray,
    target_mask: np.ndarray,
    rows: np.ndarray,
    *,
    hidden_sizes: tuple[int, int],
    epochs: int,
    batch_size: int,
    learning_rate: float,
    random_state: int,
    device: str,
    validation_rows: np.ndarray | None = None,
    patience: int = 8,
) -> tuple[object, int, float]:
    torch, nn, DataLoader, TensorDataset = _torch_modules()
    torch.manual_seed(random_state)
    if hasattr(torch.backends, "mps"):
        torch.mps.manual_seed(random_state)

    class SharedMaskedMLP(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.network = nn.Sequential(
                nn.Linear(features.shape[1], hidden_sizes[0]),
                nn.ReLU(),
                nn.Linear(hidden_sizes[0], hidden_sizes[1]),
                nn.ReLU(),
                nn.Linear(hidden_sizes[1], targets.shape[1]),
            )

        def forward(self, values):
            return self.network(values)

    model = SharedMaskedMLP().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    train_indices = np.flatnonzero(rows & target_mask.any(axis=1))
    dataset = TensorDataset(
        torch.from_numpy(features[train_indices]).float(),
        torch.from_numpy(targets[train_indices]).float(),
        torch.from_numpy(target_mask[train_indices]),
    )
    generator = torch.Generator().manual_seed(random_state)
    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        generator=generator,
        drop_last=False,
    )
    best_epoch = epochs
    best_validation = float("nan")
    best_state = None
    stale_epochs = 0
    for epoch in range(1, epochs + 1):
        model.train()
        for batch_features, batch_targets, batch_mask in loader:
            batch_features = batch_features.to(device)
            batch_targets = batch_targets.to(device)
            batch_mask = batch_mask.to(device)
            optimizer.zero_grad()
            predictions = model(batch_features)
            task_losses = []
            for task_index in range(targets.shape[1]):
                observed = batch_mask[:, task_index]
                if observed.any():
                    residuals = (
                        predictions[observed, task_index] - batch_targets[observed, task_index]
                    )
                    task_losses.append((residuals**2).mean())
            torch.stack(task_losses).mean().backward()
            optimizer.step()

        if validation_rows is None:
            continue
        model.eval()
        validation_indices = np.flatnonzero(validation_rows & target_mask.any(axis=1))
        with torch.no_grad():
            validation_features = torch.from_numpy(features[validation_indices]).float().to(device)
            validation_targets = torch.from_numpy(targets[validation_indices]).float().to(device)
            validation_mask = torch.from_numpy(target_mask[validation_indices]).to(device)
            validation_predictions = model(validation_features)
            losses = []
            for task_index in range(targets.shape[1]):
                observed = validation_mask[:, task_index]
                if observed.any():
                    losses.append(
                        (
                            (
                                validation_predictions[observed, task_index]
                                - validation_targets[observed, task_index]
                            )
                            ** 2
                        ).mean()
                    )
            validation_loss = float(torch.stack(losses).mean().cpu())
        if best_state is None or validation_loss < best_validation - 1.0e-6:
            best_validation = validation_loss
            best_epoch = epoch
            best_state = {
                name: value.detach().cpu().clone() for name, value in model.state_dict().items()
            }
            stale_epochs = 0
        else:
            stale_epochs += 1
            if stale_epochs >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, best_epoch, best_validation


def fit_masked_multitask_mlp(
    features: np.ndarray,
    targets: np.ndarray,
    train_rows: np.ndarray,
    validation_rows: np.ndarray,
    prediction_features: np.ndarray | None = None,
    *,
    hidden_sizes: tuple[int, int] = (64, 32),
    max_epochs: int = 80,
    batch_size: int = 512,
    learning_rate: float = 1.0e-3,
    patience: int = 8,
    random_state: int = 42,
    device: str = "cpu",
) -> MaskedMLPResult:
    """Tune epoch count on validation rows, then refit on all supplied rows.

    Loss is the unweighted mean of per-task MSE values after task-wise target
    standardization. Missing values never contribute to a task loss.
    """
    torch, _, _, _ = _torch_modules()
    features = np.asarray(features, dtype=np.float32)
    targets = np.asarray(targets, dtype=float)
    train_rows = np.asarray(train_rows, dtype=bool)
    validation_rows = np.asarray(validation_rows, dtype=bool)
    if targets.ndim != 2 or features.shape[0] != targets.shape[0]:
        raise ValueError("features and two-dimensional targets must share row count")
    if train_rows.shape != (len(features),) or validation_rows.shape != (len(features),):
        raise ValueError("row masks must match feature row count")
    if np.any(train_rows & validation_rows):
        raise ValueError("training and validation rows must be disjoint")
    if device != "cpu" and not (
        device == "mps" and hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
    ):
        raise ValueError(f"Unsupported or unavailable device: {device}")

    tuning_rows = train_rows | validation_rows
    tuning_means, tuning_scales = _fit_scaler(targets, train_rows)
    tuning_mask = np.isfinite(targets)
    tuning_targets = np.zeros_like(targets, dtype=np.float32)
    tuning_targets[tuning_mask] = np.broadcast_to(
        (1.0 / tuning_scales), targets.shape
    )[tuning_mask] * (targets - tuning_means)[tuning_mask]
    _, best_epoch, validation_loss = _train_epochs(
        features,
        tuning_targets,
        tuning_mask,
        train_rows,
        hidden_sizes=hidden_sizes,
        epochs=max_epochs,
        batch_size=batch_size,
        learning_rate=learning_rate,
        random_state=random_state,
        device=device,
        validation_rows=validation_rows,
        patience=patience,
    )

    final_means, final_scales = _fit_scaler(targets, tuning_rows)
    final_mask = np.isfinite(targets)
    final_targets = np.zeros_like(targets, dtype=np.float32)
    final_targets[final_mask] = np.broadcast_to(
        (1.0 / final_scales), targets.shape
    )[final_mask] * (targets - final_means)[final_mask]
    final_model, _, _ = _train_epochs(
        features,
        final_targets,
        final_mask,
        tuning_rows,
        hidden_sizes=hidden_sizes,
        epochs=best_epoch,
        batch_size=batch_size,
        learning_rate=learning_rate,
        random_state=random_state,
        device=device,
    )
    prediction_features = (
        features
        if prediction_features is None
        else np.asarray(prediction_features, dtype=np.float32)
    )
    final_model.eval()
    prediction_batches = []
    with torch.no_grad():
        for start in range(0, len(prediction_features), 8192):
            batch = torch.from_numpy(prediction_features[start : start + 8192]).float().to(device)
            prediction_batches.append(final_model(batch).cpu().numpy())
    scaled_predictions = np.concatenate(prediction_batches, axis=0)
    predictions = scaled_predictions * final_scales + final_means
    return MaskedMLPResult(
        predictions=predictions,
        best_epoch=best_epoch,
        validation_loss=validation_loss,
        training_rows=int((tuning_rows & final_mask.any(axis=1)).sum()),
        validation_rows=int((validation_rows & final_mask.any(axis=1)).sum()),
        task_observations=tuple(
            int((tuning_rows & final_mask[:, task]).sum()) for task in range(targets.shape[1])
        ),
        device=device,
    )
