# Helper to check file counts and raise appropriate errors/warnings
    def check_file_status(obj, name, sub, ses, extra_params="", is_warning_only=False):
        count = len(obj)
        if count == 0:
            msg = f"No {name} found for sub-{sub:02d}, ses-{ses:02d}{extra_params}."
            if is_warning_only:
                warnings.warn(msg)
            else:
                raise ValueError(msg)
        elif count > 1:
            msg = f"Found more than one {name}:\n{obj}.\nPlease verify your data import steps."
            raise ValueError(msg)

    # Define run string dynamically
    run_str = f", run-{runID:02d}" if use_run else ""
    task_str = f", task-{task}" if task else ""
    space_str = f", space-{space}" if space else ""
    acq_str = f", acq-{acqID}" if acqID else ""

    # Missing file checks (Raise ValueError)
    check_object(log_object, "log file", subID, sesID, f"{task_str}, {acq_str}, {run_str}")
    check_object(bold_object, "bold file", subID, sesID, f"{task_str}, {space_str}, {acq_str}, {run_str}")
    check_object(mask_object, "mask file", subID, sesID, f"{task_str}, {space_str}, {acq_str}, {run_str}")
    check_object(T1w_object, "T1w file", subID, anatID) # anatID logic usually handled in config, kept simple here
    check_object(conf_object, "confounds file", subID, sesID, f"{task_str}, {acq_str}, {run_str}")
    check_object(reg_object, "TAPAS regressors file", subID, sesID, f"{acq_str}, {run_str}", is_warning_only=True)
    check_object(movpar_object, "movement parameters file", subID, sesID, f"{task_str}, {acq_str}, {run_str}")
    check_object(out_object, "outliers file", subID, sesID, f"{task_str}, {acq_str}, {run_str}")
    check_object(T1w_to_MNI_object, "T1w to MNI transform file", subID, sesID, f"{task_str}, {space_str}, {acq_str}, {run_str}")
    
    # Transform files
    if len(orig_to_boldref_object) == 0 or len(boldref_to_T1w_object) == 0:
        raise ValueError(
            f"No orig_to_boldref to boldref_to_T1w transform files found for sub-{subID:02d}, ses-{sesID:02d}"
            f"{task_str}, {space_str}, {acq_str}, {run_str}"
        )

    # Warnings for multiple files
    if len(T1w_object) > 1:
        warnings.warn(f"Multiple anatomical files found: {[Path(to).name for to in T1w_object]}")
    
    if len(orig_to_boldref_object) > 1 or len(boldref_to_T1w_object) > 1:
        names_orig = [Path(otbo).name for otbo in orig_to_boldref_object]
        names_bold = [Path(btto).name for btto in boldref_to_T1w_object]
        warnings.warn(
            f"Multiple transformation files found: \n"
            f" * orig_to_boldref: {names_orig} \n"
            f" * boldref_to_T1w: {names_bold}")