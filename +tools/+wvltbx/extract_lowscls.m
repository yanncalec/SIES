function Sa = extract_lowscls(S, N)
    
    if nargin < 2
        Sa = S;
    else
        Sa = cell(1,length(S));
        Sa{1} = S{1};

        for n=1:length(S)-1
            if n<=N
                Sa{n+1}=S{n+1};
            else
                Sa{n+1}=zeros(size(S{n+1}));
            end
        end
    end
end