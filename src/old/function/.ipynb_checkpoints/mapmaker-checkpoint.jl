function binned_mapmaker(ψs, result_d)
    d_vector = [result_d[1], result_d[2] / 2.0, result_d[3] / 2.0]
    h2 = mean(exp.(ψs*2.0*im))
    h4 = mean(exp.(ψs*4.0*im))
    h_matrix = [1.0 0.5 * h2 0.5 * conj(h2);
        0.5*h2 0.25 * h4 0.25;
        0.5*conj(h2) 0.25 0.25 * conj(h4)
    ]
    result = h_matrix \ d_vector
    return [real(result[1]), real(result[2]), imag(result[3])]
end