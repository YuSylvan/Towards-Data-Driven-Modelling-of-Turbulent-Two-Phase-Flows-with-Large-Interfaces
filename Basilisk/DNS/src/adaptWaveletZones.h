struct AdaptZones
{
    scalar * slist;
    double * max;
    int * maxLevels;
    int * minLevels;
    int (*zone) (double, double, double);
    scalar * list;
};

astats adaptWaveletZones(struct AdaptZones p)
{
    if (p.list == NULL)
    {
        p.list = all;
    }

    if (is_constant(cm))
    {
        restriction (p.slist);
    }
    else
    {
        scalar * listr = list_concat ({cm}, p.slist);
        restriction (listr);
        free (listr);
    }

    astats st = {0, 0};

    scalar * listc = NULL;

    for (scalar s in p.list)
    {
        if (!is_constant(s) && s.restriction != no_restriction)
        {
            listc = list_add (listc, s);
        }
    }

    tree->refined.n = 0;

    static const int refined = 1 << user, too_fine = 1 << (user + 1);

    foreach_cell()
    {
        int zonei = p.zone(x, y, z);

        int cellMaxLevel = p.maxLevels[zonei];
        int cellMinLevel = p.minLevels[zonei];

        if (is_active(cell))
        {
            static const int too_coarse = 1 << (user + 2);

            if (is_leaf(cell))
            {
                if (cell.flags & too_coarse)
                {
                    cell.flags &= ~too_coarse;

                    refine_cell(point, listc, refined, &tree->refined);

                    st.nf++;
                }
                continue;
            }
            else
            {
                if (cell.flags & refined)
                {
                    // Cell has already been refined, skip its children

                    cell.flags &= ~too_coarse;

                    continue;
                }

                // Check whether the cell or any of its children is local

                bool local = is_local(cell);

                if (!local)
                {
                    foreach_child()
                    {
                        if (is_local(cell))
                        {
                            local = true;
                            break;
                        }
                    }
                }
                if (local)
                {
                    int i = 0;

                    static const int just_fine = 1 << (user + 3);

                    for (scalar s in p.slist)
                    {
                        double max = p.max[i++], sc[1 << dimension];

                        int c = 0;

                        foreach_child()
                        {
                            sc[c++] = s[];
                        }

                        s.prolongation (point, s);

                        c = 0;

                        foreach_child()
                        {
                            double e = fabs(sc[c] - s[]);

                            if (e > max && level < cellMaxLevel)
                            {
                                cell.flags &= ~too_fine;
                                cell.flags |= too_coarse;
                            }
                            else if
                            (
                                (e <= max / 1.5 || level > cellMaxLevel) &&
                               !(cell.flags & (too_coarse | just_fine))
                            )
                            {
                                if (level >= cellMinLevel)
                                {
                                    cell.flags |= too_fine;
                                }
                            }
                            else if (!(cell.flags & too_coarse))
                            {
                                cell.flags &= ~too_fine;
                                cell.flags |= just_fine;
                            }

                            s[] = sc[c++];
                        }
                    }

                    foreach_child()
                    {
                        cell.flags &= ~just_fine;

                        if (!is_leaf(cell))
                        {
                            cell.flags &= ~too_coarse;
                            if (level >= cellMaxLevel)
                            {
                                cell.flags |= too_fine;
                            }
                        }
                        else if (!is_active(cell))
                        {
                            cell.flags &= ~too_coarse;
                        }
                    }
                }
            }
        }
        else
        {
            continue;
        }
    }

    mpi_boundary_refine (listc);

    for (int l = depth(); l >= 0; l--)
    {
        foreach_cell()
        {
            int zonei = p.zone(x, y, z);

            if (!is_boundary(cell) && zonei > 0)
            {
                if (level == l)
                {
                    if (!is_leaf(cell))
                    {
                        if (cell.flags & refined)
                        {
                            // Cell was refined previously, unset the flag

                            cell.flags &= ~(refined | too_fine);
                        }
                        else if (cell.flags & too_fine)
                        {
                            if (is_local(cell) && coarsen_cell (point, listc))
                            {
                                st.nc++;
                            }

                            cell.flags &= ~too_fine;
                        }
                    }

                    if (cell.flags & too_fine)
                    {
                        cell.flags &= ~too_fine;
                    }

                    else if (level > 0 && (aparent(0).flags & too_fine))
                    {
                        aparent(0).flags &= ~too_fine;
                    }

                    continue;
                }
                else if (is_leaf(cell))
                {
                    continue;
                }
            }
        }

        mpi_boundary_coarsen(l, too_fine);
    }

    free (listc);

    mpi_all_reduce(st.nf, MPI_INT, MPI_SUM);
    mpi_all_reduce(st.nc, MPI_INT, MPI_SUM);

    if (st.nc || st.nf)
    {
        mpi_boundary_update(p.list);
    }

    return st;
}

