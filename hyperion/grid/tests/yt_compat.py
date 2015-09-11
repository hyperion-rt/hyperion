def get_frb(prj, field):
    if hasattr(prj, 'frb'):
        return prj.frb[field].value
    else:
        if prj._frb is None:
            prj._recreate_frb()
        return prj._frb[field]
