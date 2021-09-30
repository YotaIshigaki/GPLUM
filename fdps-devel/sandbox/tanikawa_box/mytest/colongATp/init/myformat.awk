{
    xpos[NR-1] = $2;
    ypos[NR-1] = $3;
    zpos[NR-1] = $4;
    n = NR;
}
END{
    omega0 = 0.3
    m = 3. * omega0 / (8. * 3.14159265359 * n);

    printf("0.0\n");
    printf("%d\n", n);
    for(i = 0; i < n; i++) {
        printf("1\n");
        printf("%+.16e\n", m);
        printf("%+.16e %+.16e %+.16e\n", xpos[i], ypos[i], zpos[i]);
        printf("%+.16e %+.16e %+.16e\n", 0., 0., 0.);
    }
}
