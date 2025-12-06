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
        hmax = p[end] - p[begin],
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
                if first
                    println("pramp - success: p=$(pcurrent)")
                else
                    println("pramp - success: p=$(pcurrent) + $(h)")
                end
            end
            pcurrent = ptrial
            if first
                first = false
            else
                h = min(h * hgrow, hmax)
            end
        catch e
            if verbose
                if first
                    println("pramp - error: p=$(pcurrent)")
                else
                    println("pramp - error: p=$(pcurrent) + $(h)")
                end
            end
            h = h * hdegrow
            if h < hmin || first
                rethrow(e)
            end
        end
    end
    return
end
