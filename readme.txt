This program takes output data from population synthesis code or from
tables with observational data from Limoges et al. 2015.
It takes R,RA,DEC,U,V,W and if (LIMOGES_CRIT_IS_USED) applies elemination 
criteria from the same article 
(if |x|>|y| & |x|>|z| then we don't take U into account)
