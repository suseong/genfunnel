function [isSol,acc] = check_sol_(pos,finalState,tsq)

isSol = 0;

if and(norm(pos - finalState(1)) < 1e-2, ~isempty(tsq))
    if and(isempty(find(diff(tsq) < -1e-3, 1)),isempty(find(tsq < -1e-3, 1)))
       isSol = 1;
    end
end

end
