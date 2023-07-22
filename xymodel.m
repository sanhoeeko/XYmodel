N=50;
%u0=rand(N,N)*(2*pi);
u0=rand(N,N)*(2*pi);
%u0=vortex(u0,26,26,-1);
%u0=vortex(u0,25,27,-1);
v0=zeros(N,N);
solver=VeloVerlet(u0,v0,0.1,0.1);
solver.start();
for t=1:1000
    solver.step();
    cla;
%     subplot(1,2,1);
    solver.plot(0.4,"charge");
%     subplot(1,2,2);
%     solver.nextEk();
%     plot(solver.eks);
    pause(0.01);
end
    