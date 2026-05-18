from chapkit import BaseConfig

from chapkit.api import (
    MLServiceBuilder, MLServiceInfo, PeriodType, ModelMetadata,
    AssessedStatus
)

from chapkit.artifact import ArtifactHierarchy
from chapkit.ml import ShellModelRunner


class MalariaSpatiotemporalConfig(BaseConfig):
    """Configuration for the malaria spatiotemporal CARBayesST model."""

    prediction_periods: int = 3
    n_sample: int = 1100
    burnin: int = 100
    thin: int = 10

info=MLServiceInfo(
    id="malaria-spatiotemporal",
    display_name="Malaria Spatiotemporal CARBayesST Model",
    period_type=PeriodType.monthly,
    version="0.0.1",
    description="Malaria Spatiotemporal lagged model",
    model_metadata=ModelMetadata(
        author="Searle et al",
        author_note="Integrating ST model with chapkit",
        author_assessed_status=AssessedStatus.orange,
        contact_email="croal008@umn.edu",
    )
)


app = (
    MLServiceBuilder(
        config_schema=MalariaSpatiotemporalConfig,
        hierarchy=ArtifactHierarchy(
            name="ml",
            level_labels={0: "ml_training_workspace", 1: "ml_prediction"},
        ),
        info=info,
        runner=ShellModelRunner(
            train_command="Rscript train.r --data {data_file}",
            predict_command=(
                "Rscript predict.r "
                "--historic {historic_file} "
                "--future {future_file} "
                "--output {output_file}"
            ),
            config_format="flat",
        ),
    )
    .build()
)
