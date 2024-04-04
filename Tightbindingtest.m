% Tight binding visualizer
% Press ctrl+C to exit

% makemovie: Create animation with param = 0...1
% pickstates: Pick out wavefunctions from the DOS plot

% Figure 3: DOS plot
% Figure 10+n: Wavefunction plots
% Figure 2: E slice
% Figure 1: Band structure

runtestcode = false;
if (runtestcode)
    plotDOS(Vlat,0);
    return;
end

makemovie = false;
sumoverrotations = false;
pickstates = true;
display3D = true;

frames = 1;
flakeconfigurations = 1;

if (makemovie)
    vid = VideoWriter('ctest.avi');
    open(vid);
    frames = 5;
end

if (sumoverrotations)
    flakeconfigurations = 10;
end

for param = [0:1/frames:1-1e-10]
    Elim = [-0.5: 0.01 :4.5];

    ginfo = struct;
    ginfo.kmeshsize = [4 4];
    ginfo.Gmeshsize = [21 21];
    ginfo.FlakeFFT = false;
    ginfo.pbc = true;
    [atdata, ginfo]=setupLattice(ginfo);
    [ginfo, U, eps, nall, indices, kindices] = kpointCalculation(atdata, ginfo, Elim, 1);
    Vlat = convertToMeshgrid(ginfo, Elim, nall);
    Vdata = getkspaceData(Vlat, ginfo, Elim);

    ginfo_flake = struct;
    ginfo_flake.kmeshsize = [8 8];
    ginfo_flake.Gmeshsize = [21 21];
    ginfo_flake.FlakeFFT = true;
    ginfo_flake.pbc = true;
    [atdata, ginfo_flake]=setupLattice(ginfo_flake);
    [ginfo_flake, U, eps, nall, indices, kindices] = kpointCalculation(atdata, ginfo_flake, Elim, 2);
    Vlat = convertToMeshgrid(ginfo_flake, Elim, nall);
    Vdata_flake = getkspaceData(Vlat, ginfo_flake, Elim);

    Vdata.V = Vdata.V + Vdata_flake.V;
    Vdata.V = sqrt(Vdata.V);

    %Vfolded = foldbands(Vdata, ginfo, Elim);
    %plotDOS(150);
    figure(3); plotDOS(Vdata, 0);
    figure(2); plotEslice(Vdata, [1 1.4], 10);
    %[xpath ypath tpath] = generatekpath([0 0; 2/3 1/3; 0.5 0; 0 0], 20);
    %plotkpath(xpath, ypath, tpath, Elim, {'Γ', 'Κ', 'M', 'Γ'});
    [xpath ypath tpath] = generatekpath(ginfo, [0 0; 0.5 0; 0.5 0.5; 0 0], 20);
    figure(1); plotkpath(Vdata, xpath, ypath, tpath, Elim, {'Γ', 'X', 'M', 'Γ'});

    [xpath ypath tpath] = generatekpath(ginfo_flake, [0 0; 0.5 0; 0.5 0.5; 0 0], 20);
    figure(4); plotkpath(Vdata, xpath, ypath, tpath, Elim, {'Γ', 'X', 'M', 'Γ'});

    if (makemovie)
        frame = getframe(gcf);
        writeVideo(vid, frame);
    end
    
    if (pickstates)
        pickWfns(atdata, ginfo, Elim, U, eps, indices, kindices);
    end
end

if (makemovie)
    close(vid);
end

function pickWfns(atdata, ginfo, Elim, U, eps, indices, kindices)
    while 1
        if (1)
            figure(3);
            waitforbuttonpress;
            pause(1);
            dcm = datacursormode(gcf);
            data = getCursorInfo(dcm);
            if (isempty(data))
                continue;
            end
            pos = data.Position;
            [tmp qindex] = min(sum((ginfo.qpts-[pos(1); pos(2)]).^2, 1));
            [tmp Eindex] = min(abs(Elim-pos(3)));
            i = indices(Eindex,qindex);
            ki = kindices(Eindex,qindex);
            if (i == 0)
                clf;
                continue;
            end
   
            nlayers = max(atdata.layer);
            for layer = [1:nlayers]
                plotWfns(atdata.X, atdata.Y, U(:,i,ki), eps(i,ki), atdata.layer, layer);
            end
        else
            figure(3);
            pause(0.5);
            [xp,yp] = ginput(1);
   
   
            [tmp x] = min(abs(qpts-xp));
            [tmp y] = min(abs(Elim-yp));
            i = indices(y,x);
            ki = kindices(y,x);
   
            if (i == 0)
                clf;
                continue;
            end
            
            nlayers = max(atdata.layer);
            for layer = [1:nlayers]
                plotWfns(atdata.X, atdata.Y, U(:,i,ki), eps(i,ki));
            end
            title("Wavefunction at q*Lat/2pi = " + num2str(qpts(x)*P/(2*pi)) + ", q*W/2pi = " + num2str(qpts(x)*gbl_W/(2*pi)))
        end
    end
end

function [ginfo, U, eps, nall, indices, kindices] = kpointCalculation(atdata, ginfo, Elim, layer)
    Epall = [0];
    nall = [];
    indices = [];
    kindices = [];
    U = [];
    eps = [];
    dE = (Elim(end)-Elim(1))/numel(Elim);

    disp(newline + "Diagonalizing Hamiltonian.");
    tic
    k_index = 1;
    for j = 1:prod(ginfo.kmeshsize)
        H_k = setupHamiltonian(atdata, ginfo, ginfo.kmesh(j,:));
        [U_k, eps_k] = eig(H_k);
        [U_k, eps_k] = sortEvecs(U_k, diag(eps_k));
        U(:,:,k_index) = U_k;
        eps = [eps eps_k];
        k_index = k_index+1;
    end
    toc

    disp(newline + "Collecting DOS contributions.");
    tic
    k_index = 1;
    for j = 1:prod(ginfo.kmeshsize)
        [nall_k, indices_k, ep_k] = getkdensity(atdata, ginfo, U(:,:,k_index), eps(:,k_index), Elim, 2*dE, layer);
        Epall = Epall + ep_k;
        nall = [nall nall_k];
        indices = [indices indices_k];
        kindices = [kindices k_index*ones(size(indices_k))];
        k_index = k_index+1;
    end
    toc
    
    [tmp, order] = sortrows(ginfo.qpts_lat', [1 2]);
    order = order';
    order = order(1,:);
    ginfo.qpts = ginfo.qpts(:,order);
    ginfo.qpts_lat = ginfo.qpts_lat(:,order);
    nall = nall(:,order);
    indices = indices(:,order);
    kindices = kindices(:,order);
end

function [atdata, ginfo] = setupGrapheneFlake()
    ginfo = struct;

    xcells = 5;
    ycells = 10;

    ginfo.R_supercell = [3*xcells 0; 0 sqrt(3)*ycells]; % Supercell periodicity
    ginfo.G = 2*pi*inv(ginfo.R_supercell)';
    
    ginfo.R_primitive = 0.5*[3 3; sqrt(3) -sqrt(3)];
    ginfo.G_primitive = (2*pi/3)*[1 1; sqrt(3) -sqrt(3)];

    atdata = struct('X', [], 'Y', [], 'Z', [], 'E0', [], 'band', [], 'layer', [], 'species', []);


    theta=pi/6;
    %theta = 0.00;
    rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    atdata = addAtoms(atdata, N=[xcells ycells], R=[3 0; 0 sqrt(3)], offset=[0 0],     Z=0, E=0, id=1, layer=1);
    atdata = addAtoms(atdata, N=[xcells ycells], R=[3 0; 0 sqrt(3)], offset=[1/3 0],   Z=0, E=0, id=1, layer=1);
    atdata = addAtoms(atdata, N=[xcells ycells], R=[3 0; 0 sqrt(3)], offset=[1/2 1/2], Z=0, E=0, id=1, layer=1);
    atdata = addAtoms(atdata, N=[xcells ycells], R=[3 0; 0 sqrt(3)], offset=[5/6 1/2], Z=0, E=0, id=1, layer=1);

    %addAtoms(atdata, N, R, offset, Z, E, id, layer)
    xcellsflake = 4;
    ycellsflake = 8;
    atdata = addAtoms(atdata, N=[xcellsflake ycellsflake], R=rot*[3 0; 0 sqrt(3)], offset=[0 0],     Z=0.5, E=1, id=2, layer=2);
    atdata = addAtoms(atdata, N=[xcellsflake ycellsflake], R=rot*[3 0; 0 sqrt(3)], offset=[1/3 0],   Z=0.5, E=1, id=2, layer=2);
    atdata = addAtoms(atdata, N=[xcellsflake ycellsflake], R=rot*[3 0; 0 sqrt(3)], offset=[1/2 1/2], Z=0.5, E=1, id=2, layer=2);
    atdata = addAtoms(atdata, N=[xcellsflake ycellsflake], R=rot*[3 0; 0 sqrt(3)], offset=[5/6 1/2], Z=0.5, E=1, id=2, layer=2);

    Hopping = [-1 -0.2; -0.2 -1];
    Hopping = Hopping+Hopping';
end

function [atdata, ginfo] = setupRectangularLattice(ginfo)
    xcells = 8;
    ycells = 8;

    xcellsflake = 5;
    ycellsflake = 5;

    theta=pi/6;
    rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    ginfo.R_primitive = [1 0; 0 1];
    ginfo.R_supercell = [xcells 0; 0 ycells]; % Supercell periodicity

    R_flake_primitive = rot*[1 0; 0 1];
    R_flake_supercell = rot*[xcellsflake 0; 0 ycellsflake];
    
    if (~ginfo.FlakeFFT)
        ginfo.G_primitive = 2*pi*inv(ginfo.R_primitive)';
        ginfo.G = 2*pi*inv(ginfo.R_supercell)';
    else
        ginfo.G_primitive = 2*pi*inv(R_flake_primitive)';
        ginfo.G = 2*pi*inv(R_flake_supercell)';
    end
    
    atdata = struct('X', [], 'Y', [], 'Z', [], 'E0', [], 'band', [], 'layer', [], 'species', [], 'connectedatoms', [], 'V', [], 'Hopping', []);


    atdata = addAtoms(atdata, N=[xcells ycells], R=[1 0; 0 1], offset=[0 0], Z=0, E=0, id=1, layer=1, pbc=ginfo.pbc);

    %addAtoms(atdata, N, R, offset, Z, E, id, layer)
    C = [(xcells-1)/2; (ycells-1)/2];
    Cflake = [(xcellsflake-1)/2; (ycellsflake-1)/2];
    dr = inv([1 0; 0 1])*(inv(rot)*C - Cflake);
    atdata = addAtoms(atdata, N=[xcellsflake ycellsflake], R=R_flake_primitive, offset=dr', Z=1.0, E=1, id=1, layer=2, pbc=ginfo.pbc);

    atdata.Hopping = -0.5;
    atdata.Hopping = atdata.Hopping+atdata.Hopping';
end

function [atdata, ginfo]=setupLattice(ginfo)
    
    [atdata, ginfo] = setupRectangularLattice(ginfo);

    V = [];
    for (i = [1:size(atdata.connectedatoms, 2)])
        Vcol = zeros(numel(atdata.E0), 1);
        Vcol(atdata.connectedatoms(1, i)) = 1/sqrt(2);
        Vcol(atdata.connectedatoms(2, i)) = -1/sqrt(2);
        V = [V Vcol];
    end
    atdata.V = V;

    [k_latX, k_latY] = meshgrid(([1:ginfo.kmeshsize(1)] - round(ginfo.kmeshsize(1)/2))*(1/ginfo.kmeshsize(1)), ...
        ([1:ginfo.kmeshsize(2)] - round(ginfo.kmeshsize(2)/2))*(1/ginfo.kmeshsize(2)));
    [GlatX, GlatY] =meshgrid(([1:ginfo.Gmeshsize(1)] - round(ginfo.Gmeshsize(1)/2)), ...
        ([1:ginfo.Gmeshsize(2)] - round(ginfo.Gmeshsize(2)/2)));

    k_lat = [k_latX(:) k_latY(:)];
    ginfo.kmesh = (ginfo.G*k_lat')';
    G_lat = [GlatX(:) GlatY(:)];
    ginfo.Glist = (ginfo.G*G_lat')';

    ginfo.qpts_lat = [(kron(k_lat', ones(1,size(G_lat,1))) + kron(ones(1,size(k_lat,1)), G_lat'))];
    %qpts = [(kron(kmesh', ones(1,size(G,1))) + kron(ones(1,size(kmesh,1)), G'))];
    ginfo.qpts = ginfo.G*ginfo.qpts_lat;

    disp("Total atoms = " + num2str(numel(atdata.E0)) + ", k-mesh size = " + num2str(prod(ginfo.kmeshsize)));
end

function Vdata=convertToMeshgrid(ginfo, Elim, nall)
    disp(newline + "Converting data to mesh grid.");
    Vdata = struct("X", [], "Y", [], "Z", [], "V", []);
    tic
    
    Xn = ginfo.kmeshsize(1)*ginfo.Gmeshsize(1);
    Yn = ginfo.kmeshsize(2)*ginfo.Gmeshsize(2);
    Zn = numel(Elim);

    Vdata.X = kron(ones(size(Elim)), ginfo.qpts_lat(1,:));
    Vdata.Y = kron(ones(size(Elim)), ginfo.qpts_lat(2,:));
    Vdata.Z = kron(Elim, ones(1,size(ginfo.qpts_lat, 2)));
    Vdata.V = reshape(nall, size(Vdata.Z));

    Vdata.X = reshape(Vdata.X, [Yn Xn,Zn]);
    Vdata.Y = reshape(Vdata.Y, [Yn Xn Zn]);
    Vdata.Z = reshape(Vdata.Z, [Yn Xn Zn]);
    Vdata.V = reshape(Vdata.V, [Zn Yn Xn]);
    Vdata.V = permute(Vdata.V, [2 3 1]);
    toc
end

function Vdatamesh=getkspaceData(Vdatalattice, ginfo, Elim)
    %Xp = kron(ones(size(Elim)), qpts(1,:));
    %Yp = kron(ones(size(Elim)), qpts(2,:));
    %Zp = kron(Elim, ones(1,size(qpts,2)));
    %Vp = reshape(nall', [1 numel(Elim)*size(qpts,2)]);
    disp(newline + "Transforming to reciprocal space coordinates.");
    tic

    Vdatamesh = struct("X", [], "Y", [], "Z", [], "V", []);

    [Vdatamesh.X Vdatamesh.Y Vdatamesh.Z] = meshgrid(2*[-pi:pi/40:pi],2*[-pi:pi/40:pi],Elim);
    
    %F = scatteredInterpolant(Xp', Yp', Zp', Vp');
    %toc
    %tic
    %V = F(X,Y,Z);

    Gpinv = inv(ginfo.G);

    Vdatamesh.V = interp3(Vdatalattice.X, Vdatalattice.Y, Vdatalattice.Z, Vdatalattice.V, Gpinv(1,1)*Vdatamesh.X + Gpinv(1,2)*Vdatamesh.Y, ...
        Gpinv(2,1)*Vdatamesh.X + Gpinv(2,2)*Vdatamesh.Y, Vdatamesh.Z);
    Vdatamesh.V(isnan(Vdatamesh.V)) = 0;
    toc
end

function Vdata_folded = foldbands(Vdata_unfolded, ginfo, Elim)
    disp(newline + "Folding bands.");
    tic
    mesh_V_f = [0];
    Vdata_folded = struct("X", [], "Y", [], "Z", [], "V", []);

    
    latticevecs = [0 0; 1 0; 1 1; 0 1; 0 -1; -1 -1; -1 0];
    shifts = (ginfo.G_primitive*latticevecs')';
    for i = 1:size(shifts,1)
        shift = shifts(i,:);
        V_tmp = interp3(Vdata_unfolded.X, Vdata_unfolded.Y, Vdata_unfolded.Z, Vdata_unfolded.V, ...
            Vdata_unfolded.X+shift(1), Vdata_unfolded.Y+shift(2), Vdata_unfolded.Z);
        V_tmp(isnan(V_tmp)) = 0;
        mesh_V_f = mesh_V_f + V_tmp;
    end
    
    [Vdata_folded.X Vdata_folded.Y Vdata_folded.Z] = meshgrid([-pi:pi/40:pi],[-pi:pi/40:pi],Elim);

    Vdata_folded.V = interp3(Vdata_unfolded.X, Vdata_unfolded.Y, Vdata_unfolded.Z, mesh_V_f, Vdata_folded.X, Vdata_folded.Y, Vdata_folded.Z);
    toc
end

function plotDOS(Vdata, cutoff)
    disp(newline + "Plotting DOS.");
    clf;
    tic
    if (cutoff == 0)
        isosurface(Vdata.X,Vdata.Y,Vdata.Z,Vdata.V);
        isocaps(Vdata.X,Vdata.Y,Vdata.Z,Vdata.V);
    else
        isosurface(Vdata.X,Vdata.Y,Vdata.Z,Vdata.V,cutoff);
        isocaps(Vdata.X,Vdata.Y,Vdata.Z,Vdata.V,cutoff);
    end
    toc
    title("Momentum resolved DOS")
    fontname("CMU Serif");
    fontsize(18, 'points');
    xlabel("k_{x}");
    ylabel("k_{y}");
    zlabel("Energy");
    camlight
    camlight(-80,-10)
    lighting gouraud
    
    %    contourf(qpts, Elim, sqrt(abs(nall)), 'LineColor', 'none');
%    title("Sqrt[Momentum resolved DOS]")
%    fontname("CMU Serif");
%    fontsize(18, 'points');
%    xlabel("kx");
%    ylabel("E");
%    xlim([-2*pi 2*pi]);
%    figure(4);
%    plot(Elim, abs(Epall));
%    fontname("CMU Serif");
%    fontsize(18, 'points');
%    xlabel("E");
%    ylabel("Relative density");
end

function plotEslice(Vdata, Elim, nslices)
    
    Eslices = [Elim(1) : (Elim(end)-Elim(1))/nslices : Elim(end)];
    h=slice(Vdata.X,Vdata.Y,Vdata.Z,Vdata.V,[],[],Eslices);
    set(h,'LineStyle','none');
    fontname("CMU Serif");
    fontsize(18, 'points');
    title("Energy slices");
    xlabel("k_{x}");
    ylabel("k_{y}");
    zlabel("Energy");
end

function [X, Y, T] = generatekpath(ginfo, nodes, resolution)
    X = [];
    Y = [];
    T = [];
    
    kcoords = (ginfo.G_primitive*nodes')';
    for i = 1:(size(kcoords,1)-1)
        A = kcoords(i,:);
        B = kcoords(i+1,:);
        t = [0:1/resolution:1];
        r = A'*(1-t) + B'*t;
        x = r(1,:);
        y = r(2,:);
        T = [T t+i-1];
        X = [X x];
        Y = [Y y];
    end
end

function plotkpath(Vdata, X, Y, T, Elim, names)

    V = interp3(Vdata.X,Vdata.Y,Vdata.Z,Vdata.V, ones(size(Elim))'*X, ones(size(Elim))'*Y, Elim'*ones(size(T)));
    contourf(ones(size(Elim))'*T, Elim'*ones(size(T)), V, 'LineStyle', 'none');
    set(gca,'XTick',[0:round(T(end))],'XTickLabel',names);
    ylabel("Energy");
    fontname("CMU Serif");
    fontsize(18, 'points');
    title("Band structure");
end

function atdata = addAtoms(atdata, params)
    arguments
        atdata struct
        params.N double
        params.R double
        params.offset double
        params.Z double
        params.E double
        params.id double
        params.layer double
        params.pbc logical
    end
    atoms = params.N(1)*params.N(2);
    [RX RY] = meshgrid([1:params.N(1)] - 1, [1:params.N(2)] - 1);
    RX = RX(:);
    RY = RY(:);
    P = params.R*[RX RY]';
    dR = params.R*params.offset';
    
    S = numel(atdata.X);

    atdata.X = [atdata.X P(1,:) + dR(1)];
    atdata.Y = [atdata.Y P(2,:) + dR(2)];
    atdata.Z = [atdata.Z params.Z*ones(1, atoms)];
    atdata.E0 = [atdata.E0 params.E*ones(1,atoms)];
    atdata.layer = [atdata.layer params.layer*ones(1,atoms)];
    atdata.species = [atdata.species params.id*ones(1,atoms)];

    if (params.pbc)
        for b = [1:params.N(1)-2]
            a1 = intersect(find(RX==b), find(RY==0))
            a2 = intersect(find(RX==b), find(RY==params.N(1)-1))
            assert(numel(a1)==1);
            assert(numel(a2)==1);
        
            atdata.connectedatoms = [atdata.connectedatoms [S+a1; S+a2]]
        end
    
        for b = [1:params.N(2)-2]
            a1 = intersect(find(RY==b), find(RX==0))
            a2 = intersect(find(RY==b), find(RX==params.N(2)-1))
            assert(numel(a1)==1);
            assert(numel(a2)==1);
        
            atdata.connectedatoms = [atdata.connectedatoms [S+a1; S+a2]]
        end
    end

    %if (pbc)
    %    atdata.V = []
end

function H = setupHamiltonian(atdata, ginfo, k)
    Totalsites = numel(atdata.E0);
    H = zeros(Totalsites, Totalsites); % Hamiltonian
    I = eye(Totalsites, Totalsites);
    V = zeros(Totalsites, 1); % Constraint projector
    Lindex = 1;
    Rindex = 1;
    %TODO

    gbl_W = atdata.X(Rindex)-atdata.X(Lindex);
    
    %V(Lindex) = -exp(i*gbl_kf*gbl_W);%-exp(-0*i*k*(X(Rindex)-X(Lindex)));
    %V(Rindex) = 1;
    %V = V/norm(V);
    
    global min_d;
    min_d = 1e99;
    for n = [1:Totalsites]
        for m = [1:Totalsites]
            for Rx = [-1 0 1]
                for Ry = [-1 0 1]
                    deltar = [atdata.X(n) atdata.Y(n) atdata.Z(n)] - [atdata.X(m) atdata.Y(m) atdata.Z(m)] - [(ginfo.R_supercell*[Rx Ry]')' 0];
                    dist = sqrt(sum(deltar.^2,'all'));
                    if (n == m && Rx == 0 && Ry == 0)
                        interactionstrength = atdata.E0(n);
                    else
                        interactionstrength = atdata.Hopping(atdata.species(n), atdata.species(m))*exp(-dist*dist*2.5)*exp(-i*k*deltar(:,1:2)');
                        min_d = min(min_d, dist);
                    end
                    H(n,m) = H(n,m) + interactionstrength;
                    %H(m,n) = H(m,n) + interactionstrength';
                end
            end
        end
    end
    if (~isempty(atdata.V))
        H = (I-atdata.V*atdata.V')*H*(I-atdata.V*atdata.V');
    end
    H = 0.5*(H+H');
end

function plotWfns(X, Y, psi, eps, layers, layer)
    global min_d;
    Totalsites = numel((X));
    xmax = max(X)+1;
    xmin = min(X)-1;
    ymax = max(Y)+1;
    ymin = min(Y)-1;
    
    xcoords = [xmin:0.05:xmax];
    ycoords = [ymin:0.05:ymax];

    [Xg, Yg] = meshgrid(xcoords, ycoords);
    
    Zg = Xg*0;
    mag = Xg*0;
    phase = Xg*NaN;

    for i = [1 : Totalsites]
        if (layers(i) == layer)
            xc = X(i);
            yc = Y(i);
            mag = mag + exp(-((Xg-xc).^2+(Yg-yc).^2)/(min_d*min_d*0.05)) * abs(psi(i));
            Zg = Zg + exp(-((Xg-xc).^2+(Yg-yc).^2)/(min_d*min_d*0.05)) * psi(i);
            phase((min_d*min_d*0.2-((Xg-xc).^2+(Yg-yc).^2)) > 0) = angle(psi(i));
        end
    end

    Zg = Zg/max(abs(psi));
    colors = complexColorize(Zg);

    figure(10+layer);
    %surf(Xg, Yg, mag, colors);
    surf(Xg, Yg, mag, colors, 'EdgeColor', 'none');
    xlabel("X");
    ylabel("Y");
    zlabel("Amplitude");
    fontname("CMU Serif");
    fontsize(18, 'points');
    view(0,90);
    %pcolor(xcoords, ycoords, colors);

    %figure(2);
    %n = getkdensity(X, Y, psi, eps, P);
    %contourf(abs(n));
    %caxis([-1/sqrt(Totalsites) 1/sqrt(Totalsites)]);
    caxis([-pi, pi]);
    colormap('hsv')
end

function [Uf, epsf] = sortEvecs(U, eps)
    [epsf, ind] = sort(eps);
    Uf = U(:, ind);
end

function [n, indices, Ep]=getkdensity(atdata, ginfo, U, E, Elim, Esmeaing, layer)
    Totalsites = numel(atdata.E0);
    m = [1:Totalsites];

    %U = U.*[zeros(12, 1); ones(13, 1); zeros(1,1)];

    V = exp(-i*ginfo.Glist*[atdata.X; atdata.Y]);
    
    n = [0];
    dnmax = [0];
    indices = zeros(numel(Elim), size(ginfo.Glist,1));
    Ep = [0];

    if (layer == -1)
        mask = ones(size(atdata.layer));
    else
        mask = (atdata.layer == layer);
    end

    if (numel(atdata.V) > 0)
        P = eye(Totalsites)-atdata.V*atdata.V';
    else
        P = eye(Totalsites);
    end

    for j = [1:Totalsites]
        nq = V*U(:,j);
        weight = mask * abs(U(:,j)).^2 * real(U(:,j)'*P*U(:,j)/(U(:,j)'*U(:,j)));
        nq = abs(nq).^2*weight;
        Esm = exp(-(Elim - E(j)).^2/(2*Esmeaing^2));
        dn = Esm'*nq';
        n = n + dn;
        indices(dn > dnmax) = j;
        dnmax = max(dn, dnmax);
        Ep = Ep + Esm';
    end
end

%small cell
%

function C = complexColorize(Z)
    C = [];
    colors = hsv;
    phase = angle(Z)/(2*pi) + 0.5;
    phase = max(phase, 0);
    phase = min(phase, 1);
    phase = phase*(size(hsv,1)-1) + 1;
    phase = round(phase);
    mag = abs(Z);
    mag = sqrt(mag)*2;
    mag = min(mag, 1);
    C(:,:,1) = reshape(colors(phase, 1), size(phase)).*mag;
    C(:,:,2) = reshape(colors(phase, 2), size(phase)).*mag;
    C(:,:,3) = reshape(colors(phase, 3), size(phase)).*mag;
end
