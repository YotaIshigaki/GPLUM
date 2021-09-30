{
    xmin = $1;
    ymin = $2;
    xmax = $4;
    ymax = $5;

    printf("%+e %+e %+e\n", xmin, ymin, 0.);
    printf("%+e %+e %+e\n", xmax, ymin, 0.);
    printf("%+e %+e %+e\n", xmax, ymax, 0.);
    printf("%+e %+e %+e\n", xmin, ymax, 0.);
    printf("\n\n");
}
