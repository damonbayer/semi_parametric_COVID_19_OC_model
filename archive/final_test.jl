sim_ids = [89, 93, 97, 101, 105, 109]

for sim_id in sim_ids
    println("sim_id: ", sim_id)
    sim_dict = subset(CSV.read("sim_table.csv", DataFrame), :sim_id => ByRow(x -> x == sim_id))[1,:]|> x -> Dict(names(x) .=> values(x))

    max_t = float(sim_dict["max_t"])
    seed = sim_dict["seed"]
    use_tests = sim_dict["use_tests"]
    use_seroprev = sim_dict["use_seroprev"]
    constant_R0 = sim_dict["constant_R0"]
    constant_alpha = sim_dict["constant_alpha"]
    constant_IFR = sim_dict["constant_IFR"]
    double_IFR_0 = sim_dict["double_IFR_0"]
    half_alpha_0 = sim_dict["half_alpha_0"]
    half_S_0 = sim_dict["half_S_0"]
    half_R0_0 = sim_dict["half_R0_0"] 

    my_model = bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, false)
    my_model_forecast = bayes_seird(data_new_deaths_forecast, data_new_cases_forecast, tests_forecast, data_seroprev_cases, seroprev_tests, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)
    my_model_forecast_missing = bayes_seird(missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)

    MAP_init = optimize_many_MAP(my_model, 100, 4, true)
    names_in_order = names(optimize(my_model, MAP(), LBFGS(linesearch = LineSearches.BackTracking())).values)[1]
    
    MAP_chains = augment_chains_with_forecast_samples(Chains(transpose(hcat(repeat([MAP_init[1]], 200)...)), names_in_order), my_model, my_model_forecast, "zeros")
    MAP_predictions = predict(my_model_forecast_missing, MAP_chains)

    CSV.write(resultsdir("MAP_predictions", savename("MAP_predictions", sim_dict, "csv")), DataFrame(MAP_predictions))
end

