function [sess] = get_lapInclude(sess)

% Get binary of indices for good laps
lapInclude = zeros(1,length(sess.ts));
for i = 1:sess.nlaps
    if isempty(find(sess.errTrials == i,1))
        lapInclude(sess.lapstt(i):sess.lapend(i)) = ones(1,diff([sess.lapstt(i) sess.lapend(i)])+1);
    end
end
sess.lapInclude = logical(lapInclude)';

end