select
y.yorf, yn.name, rr.mrna_in_pgdw, aft.avg_ratio, pr.tec1_99, pr.tec1_50, pr.tec1_20, aft.min_0, aft.min_15, aft.min_30, aft.min_45, aft.min_60, aft.min_90, aft.min_120, afc.nm_0_15, afc.nm_0_5, afc.nm_1_5, afc.nm_5_0, afc.nm_15_8, afc.nm_50_0, afc.nm_158, afc.nm_500
--y.yorf, yn.name, rr.mrna_in_pgdw, aft.avg_ratio, pr.tec1_99, pr.tec1_50, pr.tec1_20, aft.min_0, aft.min_15, aft.min_30, aft.min_45, aft.min_60, aft.min_90, aft.min_120, afc.nm_0_15, afc.nm_0_5, afc.nm_1_5, afc.nm_5_0, afc.nm_15_8, afc.nm_50_0, afc.nm_158, afc.nm_500, pr.tec1_f, pr.tec1_r, p.seq

from public.yorfs y
full join public.yorf_name yn on y.yorf = yn.yorf
full join public.ref_mrna rr on y.yorf = rr.yorf
full join public.alpha_time aft on y.yorf = aft.yorf
full join public.alpha_conc afc on y.yorf = afc.yorf
full join public.promoter_response pr on y.yorf = pr.yorf
full join public.promoters p on y.yorf = p.yorf

where
aft.yorf is not null
and pr.yorf is not null
and pr.tec1_50 > 0
and abs(aft.avg_ratio) > aft.avg_change

order by
aft.avg_ratio

