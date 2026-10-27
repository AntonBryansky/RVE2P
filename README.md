# RVE<sup>2</sup>P

<p align="left">
  <img src="RVE2P_logo.png" alt="RVE to Elastic Properties" width="500">
</p>

**From representative volume elements to effective elastic properties.**

RVE²P is an Abaqus Python script for computing the effective elastic properties of three-dimensional representative volume elements (RVEs) and repeating unit cells (RUCs) using periodic boundary conditions (PBCs).

The script, [`RVE2P.py`](RVE2P.py), automates boundary-node pairing, PBC equation creation, analysis submission, and post-processing. It calculates directional Young's moduli, Poisson's ratios, and shear moduli, with an optional effective-density calculation.

## Features

- Creates node sets for opposite faces, edges, and corners, together with reference points for loading.
- Runs three normal-loading cases and three shear-loading cases, with individual calculation switches.
- Performs an Abaqus data check before each analysis and checks job completion before post-processing.
- Extracts displacements (`U`) and reaction forces (`RF`) from the final field-output frame, using undeformed mesh dimensions.

## Requirements

Tested with **Abaqus 2024**.

- Abaqus/CAE with access to the Abaqus/Standard solver and the required licenses. Run the script in the Abaqus/CAE Python environment.
- A meshed, three-dimensional cell represented by a single part, with materials and sections assigned.
- A symmetric mesh with matching nodes on opposite faces.
- Material density definitions if effective density is requested.

## Getting started

1. Open Abaqus/CAE and create or open a model containing the meshed cell. The script defaults to model `Model-1` and part `RUC-HC`.
2. Assign material properties and sections appropriate for the elastic analysis. Check that the mesh satisfies the periodic matching requirements above.
3. Edit the input parameters near the top of `RVE2P.py` to match your model and available computing resources:

   ```python
   model_name = 'Model-1'
   part_name = 'RUC-HC'
   instance_name = 'RUC'

   num_CPUs = 16
   num_Domains = num_CPUs
   ```

   `instance_name` is the assembly instance created by the script. Set `num_CPUs` to suit your machine and license allocation.

4. Select the properties to calculate using the `calculate_*` switches described below. All switches are enabled in the current `RVE2P.py`. Set `calculate_p = False` if your materials do not include density.
5. Save a clean copy of the prepared model and choose a fresh Abaqus working directory for the analysis files.
6. In Abaqus/CAE, choose **File → Run Script** and select `RVE2P.py`. The script creates the assembly instance and PBCs, then runs the selected jobs sequentially.
7. After successful completion, read the properties in the Abaqus/CAE message area and inspect the generated `.odb` files.

Use a clean model copy for each run. The script creates named sets, reference points, boundary conditions, steps, and jobs; it also removes existing history-output requests and replaces model constraints when preparing shear loading.

## Calculations and output

The indices 1, 2, and 3 correspond to the global X, Y, and Z directions.

| Switch | Job | Calculated properties |
| --- | --- | --- |
| `calculate_E11` | `Job_E11` | Young's modulus `E11`; Poisson's ratios `nu12`, `nu13` |
| `calculate_E22` | `Job_E22` | Young's modulus `E22`; Poisson's ratios `nu21`, `nu23` |
| `calculate_E33` | `Job_E33` | Young's modulus `E33`; Poisson's ratios `nu31`, `nu32` |
| `calculate_G12` | `Job_G12` | Shear modulus `G12` in the XY plane |
| `calculate_G13` | `Job_G13` | Shear modulus `G13` in the XZ plane |
| `calculate_G23` | `Job_G23` | Shear modulus `G23` in the YZ plane |
| `calculate_p` | No solver job | Effective density: part mass divided by the full cell volume |

For `nuij`, index `i` denotes the loading direction and `j` denotes the transverse strain direction. The printed output uses labels such as `Exx`, `nu_xy`, and `Gxy`.

Solver files are written to the current Abaqus working directory, including `Job_*.odb` result databases and the associated input and diagnostic files. `RVE2P.py` prints the calculated properties; it does not automatically export a summary to CSV or a text file.

Use a consistent unit system. Elastic moduli follow the stress units of the model, and Poisson's ratios are dimensionless. Effective density follows the model's mass and length units; the printed `kg/m^3` label does not perform a unit conversion.

## Files

- `RVE2P.py` — Python script for imposing periodic boundary conditions, running analyses, and post-processing results.
- `model_RVE_core.cae` — Test Abaqus/CAE model for validation.

## Related paper

A. A. Bryansky and R. R. Badykov.

Full publication details and a DOI link will be added when available.
