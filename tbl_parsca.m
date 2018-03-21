load results/xeon/weakscaling_3d_M1e7p_N100 
[nthrs', ts(1,:)./ts] 
ans =
                         1                         1                         1
                         2         0.946711543271476         0.953263201811032
                         3         0.840453227351114         0.842464439218452
                         4         0.876075662280267         0.873650350805138
                         8         0.824219139841474         0.853226134392345
                        12         0.715853696252941         0.794486403117895
                        16         0.578721390801088         0.666275814565169
                        24         0.443053922476688          0.53493189989705

load results/xeon/strongscaling_3d_M1e7_N100
[nthrs', ts(1,:) ./ (ts.*nthrs')]           
ans =
                         1                         1                         1
                         2         0.852442331018557         0.972319549383345
                         3         0.768213353756526         0.868225151356352
                         4          0.73777359075064         0.875100984655786
                         8         0.582839341125886         0.719700329571871
                        12         0.520614223002422         0.670320228814692
                        16         0.436249009479647         0.580743875425513
                        24         0.302084921054559          0.41248049834266

HERE"S AFTER MAX PTS PRE BLOCK SET TO 1e4 INSTEAD OF 1e5:

>> fig_parallelscaling
1 threads:	 3d1 7.04 s 	 3d2 9.78 s
2 threads:	 3d1 4.07 s 	 3d2 5.5 s
3 threads:	 3d1 3.03 s 	 3d2 4.05 s
4 threads:	 3d1 2.25 s 	 3d2 2.93 s
8 threads:	 3d1 1.28 s 	 3d2 1.66 s
12 threads:	 3d1 1.04 s 	 3d2 1.24 s
16 threads:	 3d1 1.01 s 	 3d2 1.16 s
24 threads:	 3d1 0.896 s 	 3d2 1.06 s
        nthr       time(s)
            1       7.0437       9.7759
            2        4.068        5.499
            3       3.0313        4.052
            4       2.2464       2.9286
            8       1.2788        1.657
           12       1.0433       1.2362
           16        1.012       1.1631
           24      0.89633       1.0649
        nthr        par efficiencies
            1            1            1
            2      0.86574      0.88889
            3      0.77456       0.8042
            4       0.7839      0.83452
            8      0.68853      0.73748
           12      0.56264      0.65899
           16      0.43501      0.52534
           24      0.32743       0.3825
Warning: MATLAB has disabled some advanced graphics rendering features by
switching to software OpenGL. For more information, click <a
href="matlab:opengl('problems')">here</a>. 
1 threads:	 3d1 7.79 s 	 3d2 9.79 s
2 threads:	 3d1 7.71 s 	 3d2 11.4 s
3 threads:	 3d1 8.15 s 	 3d2 11.6 s
4 threads:	 3d1 7.87 s 	 3d2 10.9 s
8 threads:	 3d1 8.38 s 	 3d2 11.6 s
12 threads:	 3d1 9.45 s 	 3d2 13 s
16 threads:	 3d1 12.2 s 	 3d2 15.4 s
24 threads:	 3d1 15.6 s 	 3d2 17.8 s
        nthr       3d1 time(s)    3d2 time(s)
            1       7.7884        9.792
            2       7.7115       11.416
            3       8.1524       11.573
            4       7.8727       10.866
            8       8.3832        11.56
           12       9.4488       12.964
           16       12.217       15.367
           24       15.595         17.8
        nthr        par efficiencies
            1            1            1
            2         1.01      0.85772
            3      0.95536      0.84607
            4      0.98929      0.90116
            8      0.92905      0.84709
           12      0.82428      0.75534
           16      0.63751      0.63721
           24      0.49942      0.55012

A LITTLE BETTER
