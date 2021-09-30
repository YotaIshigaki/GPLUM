echo "0.0" > hoge
echo "24576" >> hoge
awk -f make_multi.awk xshf=-5 yshf=-3 vxshf=+2    vyshf=+0   pl008k.init >> hoge
awk -f make_multi.awk xshf=+0 yshf=+2 vxshf=+0    vyshf=-0.5 pl008k.init >> hoge
awk -f make_multi.awk xshf=+3 yshf=-4 vxshf=-0.75 vyshf=+1   pl008k.init >> hoge
mv hoge pl008kx3.init
