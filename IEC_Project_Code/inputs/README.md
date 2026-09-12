# Input interfaces

[schemas.R](schemas.R) defines named, zero-row data-frame objects for the tabular inputs. They document field names and basic types without including observations or synthetic results. Sourcing the file creates an `input_schemas` list and writes nothing.

Populate the appropriate tables outside the source repository and supply their location with `METHODS_INPUT_DIR`, or place them here for a local analysis. Tabular inputs are tab-separated with a header row. The repository ignores input data by default.

The schemas describe the main interfaces; optional aliases and additional analysis-specific fields are documented in the scripts. Biological group labels, missing-value conventions, units, and coordinate builds must agree with the chosen analysis. For column semantics, see [input contracts](../docs/input_contracts.md).

Prepared Seurat objects, pathway collections, BAMs, clone tables, FASTQs, genome references, and annotations are external inputs. None are bundled.

The targeted vector-assay interfaces are defined separately in [vector_assay/input_schemas.py](../vector_assay/input_schemas.py); their default data location is `inputs/vector_assay/`.
