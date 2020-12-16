BEGIN{
    n = 0;
}
{
    if(NR <= 2) {
        next;
    } else if (NR % 4 == 0) {
        m[n] = $1;
    } else if (NR % 4 == 1) {
        px[n] = $1;
        py[n] = $2;
        pz[n] = $3;
    } else if (NR % 4 == 2) {
        vx[n] = $1;
        vy[n] = $2;
        vz[n] = $3;
        n++;
    } else {
        next;
    }
}
END{
    printf("0.0\n");
    printf("%d\n", n);
    printf("3\n");
    for(i = 0; i < n; i++)
        printf("%+.16e\n", m[i]);
    for(i = 0; i < n; i++)
        printf("%+.16e %+.16e %+.16e\n", px[i], py[i], pz[i]);
    for(i = 0; i < n; i++)
        printf("%+.16e %+.16e %+.16e\n", vx[i], vy[i], vz[i]);
}

