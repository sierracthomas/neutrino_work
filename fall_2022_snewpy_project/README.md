# Instructions:
- clone or copy (at least) `scripts/` locally
- be sure to have [SNOwGLoBES](https://github.com/SNOwGLoBES/snowglobes) and [SNEWPY](https://snewpy.readthedocs.io/en/stable/gettingstarted.html#installation) installed
- start in `scripts/` folder to start by running SNEWPY
- open `generate.ipynb` and change some variables:
  - `SNOwGLoBES_path` to the path you've installed snowglobes in
  - `SNEWPY_models_base` to the path where your SNEWPY models live
  - update parameters like `detector`, `distance`, `transformation`....
  - change times to the times you'd like to run (NOTE: this varies by model, will get an error if your time does not fall within the parameters of the model)
  - in the second code-block, change the `modeltype` and `model` variables like `modeltype = matrix[x][0]` and `model = matrix[x][1]` where `x` corresponds to the row of model you want. Change the matrix if you don't see yours