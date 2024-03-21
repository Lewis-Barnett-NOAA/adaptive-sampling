--sea ice query for Lewis Barnett
--Monthly average for March in the Bering Sea
select round(avg(sea_ice_fraction),2) mean_march_ice, extract(year from read_date) year
from afsc.erddap_crw_sst a
left join afsc.erddap_crw_sst_spatial_lookup b on a.crw_id=b.id
where ecosystem='Eastern Bering Sea'
and extract(month from read_date)=3
group by extract(year from read_date)
order by extract(year from read_date);