    if (argc < 6 && master)
    {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s\n", argv[0]);
        fprintf(stderr, "    [min refinement level]\n");
        fprintf(stderr, "    [max refinement level]\n");
        fprintf(stderr, "    [refinement tolerance]\n");
        fprintf(stderr, "    [write time]\n");
        fprintf(stderr, "    [end time]\n\n");
        fflush(stderr);
        exit(0);
    }

    // Read min and max refinement levels

    sscanf(argv[1], "%d", &minLevel);
    sscanf(argv[2], "%d", &maxLevel);
    sscanf(argv[3], "%lf", &epsilon);
    sscanf(argv[4], "%lf", &writeTime);
    sscanf(argv[5], "%lf", &endTime);

    if (master)
    {
        fprintf(stdout, "Refinement level min, max = %i, %i\n", minLevel, maxLevel);
        fprintf(stdout, "Refinment epsilon = %.8g\n", epsilon);
        fprintf(stdout, "Write time = %.8g, end time = %.8g\n", writeTime, endTime);
        fflush(stdout);
    }
