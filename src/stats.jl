module Stats

using LinearAlgebra

function chi_sq(obs, pred, err)
    safe_err = map(e -> e == 0.0 ? 1.0 : e, err)
    if any(err .== 0.0)
        @warn "chi_sq: $(count(err .== 0.0)) obs have err==0; substituting err=1 to avoid divide-by-zero"
    end
    sum(((obs .- pred) .^ 2) ./ safe_err .^ 2)
end


function reduced_chi_sq(obs, pred, err, n_param)
    dof = length(obs) - n_param
    if dof <= 0
        @warn "reduced_chi_sq: n_obs ($(length(obs))) <= n_params ($n_param); returning NaN"
        return NaN
    end
    chi_sq(obs, pred, err) / dof
end


function vector_chi_sq_misfit(obs1, obs2, pred1, pred2, err1, err2, cov)
    resid1 = obs1 - pred1
    resid2 = obs2 - pred2

    cc = [err1^2 cov; cov err2^2]

    p_ =[resid1 resid2] * pinv(cc) * [resid1; resid2]
    return p_[1]
end


function reduced_chi_sq_mccaffrey(misfits, n_params; n_obs=length(misfits))
    # McCaffrey's "normalized chi" (DefNode), chi_n = sqrt(chi^2 / dof).
    # n_obs must count scalar independent observations; for entries that are
    # themselves multi-component Mahalanobis chi^2 values (e.g. fused E+N GNSS),
    # pass the true scalar obs count, not length(misfits).
    dof = n_obs - n_params

    if dof <= 0
        @warn "reduced_chi_sq_mccaffrey: n_obs ($n_obs) <= n_params ($n_params); returning NaN"
        return NaN
    end

    return sqrt(sum(misfits) / dof)
end


function akaike_info_crit()
end

function gauss_likelihood_scalar(obs, pred, err)

    term1 = 1 / (sqrt(2 * pi) * err)
    term2 = exp(- (obs - pred)^2 / (2 * err^2))

    term1 * term2
end


function gauss_likelihood_vector(obs, pred, err)

    term1 = 1 ./ (sqrt(2 * pi) .* err)
    term2 = exp.(- (obs .- pred).^2 ./ (2 .* err.^2))

    term1 .* term2
end

"""
    residual_standard_error(obs, pred, num_vars, num_params)

Unbiased estimate of residual standard deviation (dof-adjusted, not classical
RMSE). Historically also exposed as `RMSE(obs, pred, num_vars, num_params)` for
backward compatibility.
"""
function residual_standard_error(obs, pred, num_vars, num_params)
    dof = num_vars - num_params
    if dof <= 0
        @warn "residual_standard_error: num_vars ($num_vars) <= num_params ($num_params); returning NaN"
        return NaN
    end
    resids = pred .- obs
    MSE = 1 / dof * resids' * resids
    sqrt(MSE)
end

# Back-compat: old 4-arg name was RMSE even though it's really the residual
# standard error. Keep the alias but route through the correctly-named impl.
RMSE(obs, pred, num_vars, num_params) = residual_standard_error(obs, pred, num_vars, num_params)


function RMSE(obs, pred)
    resids = pred .- obs
    MSE = sum(resids.^2) / length(resids)
    sqrt(MSE)
end

function vector_misfit(obs1, obs2, pred1, pred2)
    resid1 = obs1 - pred1
    resid2 = obs2 - pred2

    return sqrt(resid1^2 + resid2^2)
end

end # module
