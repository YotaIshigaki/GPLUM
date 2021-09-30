#argv: xshf, yshf, vxshf, vyshf
{
    if(NR <= 2) {
        next;
    } else if(NR % 4 == 0) {
        printf("%+.16e\n", $1);
    } else if(NR % 4 == 1) {
        $1 += xshf;
        $2 += yshf;
        printf("%+.16e %+.16e %+.16e\n", $1, $2, $3);
    } else if(NR % 4 == 2) {
        $1 += vxshf;
        $2 += vyshf;
        printf("%+.16e %+.16e %+.16e\n", $1, $2, $3);        
    } else {
        print $0;
    }
    
}
