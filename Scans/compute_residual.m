function RES = compute_residual(data,sims,n_steps)

if isrow(data); data = data'; end;

datas = repmat(data,[1 n_steps]);

diff = datas - sims;

sq = diff.^2;

RES = squeeze(sum(sq));

