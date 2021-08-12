function artifacts = normalize_artifacts(artifacts)
% Takes a 2-column array of artifacts, sorts it by start time, and modifies it if necessary
% so that there are no overlaps without changing the segments of data that are marked as artifact.

if isempty(artifacts)
    artifacts = zeros(0, 2);
    return;
end

assert(all(diff(artifacts, [], 2) > 0), 'Invalid artifacts (of negative or zero length)');
artifacts = sortrows(artifacts);

kA = 2;
while kA <= size(artifacts, 1)
    if artifacts(kA, 1) > artifacts(kA-1, 2)
        % Normal case, there's a good segment between this artifact and the previous.
        kA = kA + 1;
    elseif artifacts(kA, 2) <= artifacts(kA-1, 2)
        % In this case, this artifact is completely encompassed by the previous one. Delete.
        artifacts(kA, :) = [];
    else
        % This artifact begins within the previous one but extends beyond it. Combine them.
        artifacts(kA-1, 2) = artifacts(kA, 2);
        artifacts(kA, :) = [];
    end
end

end
