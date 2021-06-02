###### x: $3, y: $4
{
    if(NR == 1) {
        nx = int((xmax - xmin) / width);
        ny = int((ymax - ymin) / width);
        for(i = 0; i < nx; i++) {
            for(j = 0; j < ny; j++) 
                den[i,j] = 0.;
        }
    }

    if((xmin <= $3 && $3 < xmax) && (ymin <= $4 && $4 < ymax)) {
        ix = int(($3 - xmin) / width);
        iy = int(($4 - ymin) / width);
        den[ix,iy] += 1;
    }
}
END{
    ldenmax = log(den[0,0] + 1.) / log(10.);
    ldenmin = log(den[0,0] + 1.) / log(10.);
    for(i = 0; i < nx; i++) {
        for(j = 0; j < ny; j++) {
            lden[i,j] = log(den[i,j] + 1.) / log(10.);
            if(lden[i,j] > ldenmax)
                ldenmax = lden[i,j];
            if(lden[i,j] < ldenmin)
                ldenmin = lden[i,j];
        }
    }

    ldenmaxinv = 1. / ldenmax;
    for(i = 0; i < nx; i++) {
        for(j = 0; j < ny; j++) {
            lden[i,j] *= ldenmaxinv;
        }
    }    

    for(i = 0; i < nx; i++) {
        for(j = 0; j < ny; j++) {
            x = i * width + xmin;
            y = j * width + ymin;
            printf("%+e %+e %+e\n", x, y, lden[i,j]);
        }
        printf("\n");
    }
}
