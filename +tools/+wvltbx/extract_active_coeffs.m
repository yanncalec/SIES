function S = extract_active_coeffs(W, maskc);

    S=[];

    for n = 1:length(W)
        if n==1
            S{n}=W{n};
        else
            S{n}=W{n}.*maskc{n};
        end
    end
end