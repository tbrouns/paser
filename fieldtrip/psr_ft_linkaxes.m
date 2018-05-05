function psr_ft_linkaxes(figHandle,axs,titles,climType)

if (nargin < 4); climType = []; end

set(0,'CurrentFigure',figHandle);

nAxes   = size(axs,1);
climits = zeros(nAxes,2);

for i = 1:nAxes; climits(i,:) = axs(i,1).CLim; end

climits = mean(climits);
climits = max(abs(climits)); 

if (strcmp(climType,'both')); climits = [-climits climits];
else,                         climits =        [0 climits];
end

for i = 1:nAxes
    set(figHandle,'currentaxes',axs(i,1));
    caxis(climits); 
    set(axs(i,2),'CLim',climits);
    if (length(titles) >= i); title(titles{i}); end
end

axesAll = findall(figHandle,'type','axes');
n = length(axesAll) / nAxes;
if (mod(n,1) == 0)
    for i = 1:n
        for j = 1:nAxes
            k = i + n * (j - 1);
            axesTemp(j) = axesAll(k); 
        end
        linkaxes(axesTemp);
    end
end
                
end