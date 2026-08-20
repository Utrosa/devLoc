def addNuisance(bunch, confounds_path, confounds_names):
    '''
    Parameters:
        bunch: a list with a Bunch object, created by parsing logfiles of the experimental task.
        confounds_path: The path to filtered confounds (physiological regressors files - TAPAS/BIOPAC).
        confounds_names: The column names of these confounds/regressors.

    Returns:
        Bunch: a list with a bunch object that includes the regressors.
    '''
    import pandas as pd
    from nipype.interfaces.base import Bunch

    # Select the bunch object
    design_bunch = bunch[0]

    # Read the confounds file
    all_confounds = pd.read_csv(
        confounds_path,
        sep='\t',
        header=None,
        names=confounds_names,
        index_col=False
    )

    # Remove the last row of the dataframe which correspnds to the noise scan volume
    # This is not removed for physiological data ....
    all_confounds = all_confounds.iloc[:-1]

    # Convert to the required format for SPM Bunch
    regressors = [all_confounds[col].tolist() for col in all_confounds.columns]
    
    # Add regressors
    design_bunch.regressors = regressors
    design_bunch.regressor_names = confounds_names
    
    return [design_bunch]