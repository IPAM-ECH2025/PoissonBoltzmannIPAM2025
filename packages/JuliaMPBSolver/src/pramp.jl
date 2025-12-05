"""
     pramp(
        f;
        p = [0, 1],
        hmin = (p[end] - p[begin]) / 1000,
        hmax = (p[end] - p[begin]) / 10,
        h = hmax,
        hgrow = 1.2,
        hdegrow = 0.5,
        verbose = false
    )

Run 'f(p)'  successively with parameter values from the range given by 'p'.
Stepsize is given by `h`. If `f` throws an error, solution is retried with a lower
value of `h`, otherwise, `h` is increased unless it exceeded `hmax`.
"""
function pramp(
        f;
        p = [0, 1],
        hmin = (p[end] - p[begin]) / 1000,
        hmax = (p[end] - p[begin]) / 10,
        h = hmax,
        hgrow = 1.2,
        hdegrow = 0.5,
        verbose = false
    )
    pcurrent = p[begin]
    first = true
    h = min(h, hmax)
    while pcurrent < p[end]
        try
            if first
                ptrial = pcurrent
            else
                ptrial = min(pcurrent + h, p[end])
            end
            f(ptrial)
            if verbose
                println("pramp - success: p=$(pcurrent), h=$(h)")
            end
            pcurrent = ptrial
            h = min(h * hgrow, hmax)
            first = false
        catch e
            if verbose
                println("pramp - error: p=$(pcurrent), h=$(h)")
            end
            h = h * hdegrow
            if h < hmin || first
                rethrow(e)
            end
        end
    end
    return
end
