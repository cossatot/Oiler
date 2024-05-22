module Stats

function chi_sq(obs, pred, err)
    sum(((obs .- pred).^2) ./ err.^2)
end


function reduced_chi_sq(obs, pred, err, n_param)
    chi_sq(obs, pred, err) / (length(obs) - n_param)
end


function vector_misfit(obs1, obs2, pred1, pred2, err1, err2, cov)
    resid1 = obs1 - pred1
    resid2 = obs2 - pred2

    cc = [err1^2 cov; cov err2^2]

    p_ =[resid1 resid2] * cc * [resid1; resid2]
    return p_[1]
end


function reduced_chi_sq_mccaffrey(misfits, n_params)
    n_obs = length(misfits)

    return sqrt(sum(misfits) / (n_obs - n_params))

end


function akaiki_info_crit()
end

function gauss_likelihood_scalar(obs, pred, err)

    term1 = 1 / sqrt(2 * pi) * err
    term2 = exp(- (obs - pred)^2 / (2 * err^2))
    
    term1 * term2
end


function gauss_likelihood_vector(obs, pred, err)

    term1 = 1 / sqrt(2 * pi) .* err
    term2 = exp.(- (obs .- pred).^2 / (2 * err.^2))
    
    term1 .* term2
end

function RMSE(obs, pred, num_vars, num_params)
    resids = pred .- obs
    MSE = 1 / (num_vars - num_params) * resids' * resids
    sqrt(MSE)
end


end # module