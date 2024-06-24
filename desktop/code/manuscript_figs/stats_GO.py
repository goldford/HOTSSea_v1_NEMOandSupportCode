import numpy as np


# alt pearsons R
def pearson2(obs, mod):
    mod = np.asarray(mod)
    obs = np.asarray(obs)
    obs = np.asarray(obs)
    R = (1 / (np.nanstd(mod) * np.nanstd(obs))) * (np.sum((mod - np.nanmean(mod)) * (obs - np.nanmean(obs))) / len(obs))
    return R


def willmott1981(obs, mod, axis=0):
    mod = np.asarray(mod)  # add by GO
    obs = np.asarray(obs)

    num = np.nansum((mod - obs) ** 2, axis=axis)
    obs_mean = np.nanmean(obs, axis=axis)
    dM = np.abs(mod - obs_mean)
    dO = np.abs(obs - obs_mean)
    den = np.nansum((dM + dO) ** 2, axis=axis)
    if den == 0:
        return np.nan
    else:
        return np.max([0, 1 - num / den])

