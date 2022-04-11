if !@isdefined(max_t)
    error("Must define max_t before loading data")
end

if !@isdefined(max_t_forecast)
    @warn("max_t_forecast is not defined. Assigning max_t_forecast = max_t")
    max_t_forecast = max_t 
end

full_dat = CSV.read("data/oc_data.csv", DataFrame)[:, [:time, :cases, :tests, :deaths]]
seroprev_dat = CSV.read("data/oc_seroprev_data.csv", DataFrame)

max_index = searchsortedlast(full_dat[!, :time], max_t)
max_index_forecast = searchsortedlast(full_dat[!, :time], max_t_forecast)

max_index_seroprev = searchsortedlast(seroprev_dat[!, :time], max_t)
max_index_seroprev_forecast = searchsortedlast(seroprev_dat[!, :time], max_t_forecast)

obstimes = float.(full_dat[1:max_index, :time])
obstimes_forecast = float.(full_dat[1:max_index_forecast, :time])

seroprev_times = seroprev_dat[1:max_index_seroprev, :time]
seroprev_times_forecast = seroprev_dat[1:max_index_seroprev_forecast, :time]

param_change_times = obstimes[1:(end - 1)]
param_change_times_forecast = obstimes_forecast[1:(end - 1)]

tests = full_dat[1:max_index, :tests]
tests_forecast = full_dat[1:max_index_forecast, :tests]

seroprev_tests = seroprev_dat[1:max_index_seroprev, :seroprev_tests]
seroprev_tests_forecast = seroprev_dat[1:max_index_seroprev_forecast, :seroprev_tests]

data_new_deaths = full_dat[1:max_index, :deaths]
data_new_deaths_forecast = full_dat[1:max_index_forecast, :deaths]
missing_new_deaths_forecast = repeat([missing], max_index_forecast)

data_new_cases = full_dat[1:max_index, :cases]
data_new_cases_forecast = full_dat[1:max_index_forecast, :cases]
missing_new_cases_forecast = repeat([missing], max_index_forecast)

data_seroprev_cases = seroprev_dat[1:max_index_seroprev, :seroprev_cases]
data_seroprev_cases_forecast = seroprev_dat[1:max_index_seroprev_forecast, :seroprev_cases]
missing_seroprev_cases_forecast = repeat([missing], max_index_seroprev_forecast)