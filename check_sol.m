function [isSol,acc] = check_sol(pos,initState,finalState,tsq,actLimit)

isSol = 0;
acc = [];

if and(norm(pos - finalState(1)) < 1e-2, ~isempty(tsq))
    if and(isempty(find(diff(tsq) < -1e-3, 1)),isempty(find(tsq < -1e-3, 1)))
        acc = calc_acc(initState(3),tsq,actLimit(1));
        if max(abs(acc)) < abs(actLimit(2)) + 1e-6
            isSol = 1;
        end
    end
end

end
