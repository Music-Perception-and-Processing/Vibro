function qf_bartlett(x)
% Bartlett's test of sphericity 
    R = corrcoef(x);  % R:  correlation matrix
    n = size(x,1);  % n:  sample size
    detR = det(R);
    p = size(R,1);
    chisq_statistic = -log(detR) * (n - 1 - (2 * p + 5)/6);
    df = p * (p - 1)/2;  % degrees of freedom
    pval = 1-chi2cdf(chisq_statistic,df);  
    res ="=============================="+newline+...
    "Bartlett's test of sphericity" + newline+...
    "Chi-square statistic: "+chisq_statistic+newline+...
        "Degrees of freedom: " + df+newline+...
        "P-value: "+pval+newline+...
        "=============================="+newline;
    fprintf(res)
end