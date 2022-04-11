function seirdc_log_ode!(du, u, p, t)
    (S, E, I, R, D, C) = exp.(u)
    (β, γ, ν, IFR) = p
    N = S + E + I + R + D
  
    infection = β * I * S / N
    progression = γ * E
    recovery = ν * (1 - IFR) * I
    death = ν * IFR * I
  
    @inbounds begin
        du[1] = -infection / S # S
        du[2] = (infection - progression) / E # E
        du[3] = (progression - (recovery + death)) / I # I
        du[4] = recovery / R # R
        du[5] = death / D # D
        du[6] = progression / C # Cumulative Progressions
    end
    nothing
end