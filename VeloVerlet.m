classdef VeloVerlet < handle
    properties
        ji
        dt
        u
        v
        m
        q
        eks
        qs
    end
    methods
        function self=VeloVerlet(u0,v0,ji,dt)

% input:    u: initial angle phi
%           v: initial d/dt phi
%           ji:=J/2I with
%               -> J: coupling constant (positive for ferromagnetism)
%               -> I: inertia
%           dt: time step

            self.ji=ji;
            self.dt=dt;
            self.u=u0;
            self.v=v0;
            self.eks=[];
        end

        function step(self)

% input:    --from self
% output:   U: phi(t+dt)
%           V: \dot{phi}(t+dt)

            U=self.u+self.v*self.dt+self.m*(self.dt*self.dt);
            M=calForce(self.ji,U);
            V=self.v+(self.m+M)*self.dt;
% update cache
            self.u=U;self.v=V;self.m=M;
        end

        function start(self)
            self.m=calForce(self.ji,self.u);
        end

        function [x,y]=toVectorField(self)
            x=cos(self.u);
            y=sin(self.u);
        end

        function plot(self,scale_,scattertype)
            [x,y]=self.toVectorField();
            [n1,n2]=size(self.u);
            [X0,Y0]=meshgrid(1:n1,1:n2);

% slightly shift the start points in order to align the centers of grids.
            X=X0-x.*scale_/2; Y=Y0-y.*scale_/2;

            hold on
            if scattertype=="angle"
                scatter(reshape(X0,[n1*n2,1]), ...
                        reshape(Y0,[n1*n2,1]), ...
                        100*scale_, ...
                        phaseColor(self.u), ...
                        "filled");

            elseif scattertype=="charge8"
                Q=self.getTopoChargeField();
                scatter(reshape(X0,[n1*n2,1]), ...
                        reshape(Y0,[n1*n2,1]), ...
                        100*scale_, ...
                        reshape(Q,[n1*n2,1]), ...
                        "filled");
                clim([-2,2]);
                colormap jet
                colorbar;

            elseif scattertype=="charge"
                Q=self.getTopoChargeFieldSimple();
                scatter(reshape(X0,[n1*n2,1]), ...
                        reshape(Y0,[n1*n2,1]), ...
                        100*scale_, ...
                        reshape(Q,[n1*n2,1]), ...
                        "filled");
                clim([-2,2]);
                colormap jet
                colorbar;
            end

            quiver(X,Y,x,y,scale_,"black");
            axis equal
        end

        function Ek=getKineticEnergy(self)
            Ek=sum(self.v.*self.v,"all");
        end

        function nextEk(self)
            self.eks=[self.eks,self.getKineticEnergy()];
        end

        function Q=getTopoChargeField(self)
            u_main=mod(self.u/(2*pi),1)*(2*pi);
            e=circshift(u_main,1,1);
            w=circshift(u_main,-1,1);
            n=circshift(u_main,1,2);
            s=circshift(u_main,-1,2);
            ne=circshift(n,1,1);
            nw=circshift(n,-1,1);
            se=circshift(s,1,1);
            sw=circshift(s,-1,1);
            dphi1=minAbsPhase(ne-e);
            dphi2=minAbsPhase(n-ne);
            dphi3=minAbsPhase(nw-n);
            dphi4=minAbsPhase(w-nw);
            dphi5=minAbsPhase(sw-w);
            dphi6=minAbsPhase(s-sw);
            dphi7=minAbsPhase(se-s);
            dphi8=minAbsPhase(e-se);
            Q=dphi1+dphi2+dphi3+dphi4+dphi5+dphi6+dphi7+dphi8;
        % cache
            self.q=Q;
        end

        function Q=getTopoChargeFieldSimple(self)
            u_main=mod(self.u/(2*pi),1)*(2*pi);
            e=circshift(u_main,1,1);
            w=circshift(u_main,-1,1);
            n=circshift(u_main,1,2);
            s=circshift(u_main,-1,2);
            dphi1=minAbsPhase(e-n);
            dphi2=minAbsPhase(n-w);
            dphi3=minAbsPhase(w-s);
            dphi4=minAbsPhase(s-e);
            Q=dphi1+dphi2+dphi3+dphi4;
        % cache
            self.q=Q;
        end

        function rho=getDefectDensity(self,sign)
            if sign==1
                Qsum=sum(self.q(self.q>0));
            elseif sign==-1
                Qsum=sum(self.q(self.q<0));
            end
            [n1,n2]=size(self.q);
            rho=Qsum/(n1*n2);
        end

        function nextQ(self)
            self.qs=[self.qs, ...
                [self.getDefectDensity(1);
                 self.getDefectDensity(-1)]
            ];
        end

    end
end

function M=calForce(ji,U)
    M=-ji*(sin(U-circshift(U,1,1))+...
           sin(U-circshift(U,-1,1))+...
           sin(U-circshift(U,1,2))+...
           sin(U-circshift(U,-1,2)) ...
           );
end

function rgb=phaseColor(phi)
    [n1,n2]=size(phi);
    h=reshape(phi/(2*pi),[n1*n2,1]);
    h=mod(h,1);
    s=ones(size(h))*0.5;
    v=ones(size(h));
    hsv=[h,s,v];
    rgb=hsv2rgb(hsv);
end

function mu=minAbsPhase(u) % restrict the phase difference in [-pi,pi]
    mu=mod(u/(2*pi)+0.5,1)-0.5;
end