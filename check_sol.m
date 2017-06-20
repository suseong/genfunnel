function isSol = check_sol(pos,initState,finalState,tsq,actLimit)

isSol = 0;
if and(norm(pos - finalState(1)) < 1e-2, ~isempty(tsq))
    if isempty(find(diff(tsq) < 0, 1))
        acc = calc_acc(initState(3),tsq,actLimit(1));
        if max(abs(acc)) < abs(actLimit(2)) + 1e-1
            isSol = 1;
        end
    end
end

end
