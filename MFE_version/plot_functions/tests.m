%a bunch of tests to make sure all of the functions are working properly

parent_dir = cd(cd('..'));
addpath(parent_dir);

%set the dimension for all the tests
dim = 10;

%let's make a lot of bodies that we'll test
[P_cube, cube_vol] = makeBody('cube',dim);
[P_long_box, long_box_vol] = makeBody('long_box',dim);
[E_ellipsoid, ellipsoid_vol] = makeBody('ellipsoid',dim);
[E_ball, ball_vol] = makeBody('ball',dim);
[P_simplex_iso, simplex_iso_vol] = makeBody('isotropic_simplex',dim);
[P_simplex_std, simplex_std_vol] = makeBody('standard_simplex',dim);
[P_birkhoff] = makeBody('birkhoff', ceil(sqrt(dim)));
[P_triply_stoch] = makeBody('triply_stochastic', ceil(dim^(1/3)));
P_halfspace.A = [-1 zeros(1,dim-1)];
P_halfspace.b = 0;
P_halfspace.p = [0.5; zeros(dim-1,1)];
P_long_box2 = P_long_box;
P_long_box2.p = [-90; .5*ones(dim-1,1)];

commands = [    
    '[v]=Volume(P_cube);                                                             ';
    '[v]=Volume(P_long_box,[],.2,''-round'');                                          ';
    '[v]=Volume(P_simplex_iso);                                                      ';
    '[v]=Volume(P_simplex_std,[],.2,''-round'');                                       ';
    '[v]=Volume([],E_ball);                                                          ';
    '[v]=Volume([],E_ellipsoid,.2,''-round'');                                         ';
    '[v]=Volume(P_halfspace,E_ball);                                                 ';
    '[v]=Volume(P_birkhoff);                                                         ';
    '[v]=Volume(P_triply_stoch);                                                     ';
    '[v]=Volume(P_long_box2,[],.2,''-round'');                                         ';
    '[v]=Volume(P_cube,[],.2,''-walk har'');                                           ';
    '[v]=Volume(P_cube,[],.2,''-walk ball'');                                          ';
    '[v]=Volume(P_cube,[],.2,''-walk har -C 8 -verb 2 -min_st 1000 -round -a_stop 0'');';
    '[v]=Volume(P_cube,[],.2,''-plot'');                                               ';
    '[v]=Volume(P_cube,[],.2,''-plot_c'');                                             ';
    '[v]=Cube(dim,.2,1);                                                             ';
    '[v]=Cube(dim,.2,2);                                                             ';
    '[v]=Cube(dim,.2,3);                                                             ';
    '[v]=Simplex(dim,.2,1);                                                          ';
    '[v]=Simplex(dim,.2,2);                                                          ';
    'Sample(P_cube);                                                                 ';
    'Sample(P_simplex_iso);                                                          ';
    'Sample(P_simplex_std,[],.2,''-round'');                                           ';
    'Sample([],E_ball);                                                              ';
    'Sample([],E_ellipsoid,.2,''-round'');                                             ';
    'Sample(P_halfspace,E_ball);                                                     ';
    
    ];

outputs = [...
    cube_vol;
    long_box_vol;
    simplex_iso_vol;
    simplex_std_vol;
    ball_vol;
    ellipsoid_vol;
    ball_vol/2;
    -1;
    -1;
    long_box_vol;
    cube_vol;
    cube_vol;
    cube_vol;
    cube_vol;
    cube_vol;
    cube_vol;
    -1;
    cube_vol;
    simplex_iso_vol;
    simplex_std_vol;
    -1;
    -1;
    -1;
    -1;
    -1;
    -1
    ];

errors = zeros(size(commands,1),1);

for i=1:size(commands,1)
   eval(commands(i,:));
   if outputs(i)~=-1
       errors(i) = abs(outputs(i)-v)/outputs(i);
      if errors(i) > 0.8
         error('The error on the command %s is suspiciously high (%f).\n', strtrim(commands(i,:)), errors(i)); 
      end
   end
end

for i=1:size(commands,1)
   fprintf('Test case %d: %s', i, strtrim(commands(i,:)));
   if outputs(i)==-1
      fprintf('...no target answer, but completed successfully.\n');
   else
       fprintf('...completed with eps=%f.\n', errors(i));
   end
end

fprintf('Highest error was %f\n', max(errors));