
%%%%%%%%===================================================================
%%%%%%%Typhoon_data_new:1.Typhoon_rank
%%%%%%%                2.Typhoon occur time (days since 1858-11-17 00:00:00,matlab_time)
%%%%%%%                3.the longitude of Typhoon center
%%%%%%%                4.the latitude of Typhoon center
%%%%%%%                5.maximum sustained wind velocity(m/s,has been transferred)
%%%%%%%                6.Radius of maximum winds (not best tracked,nmile)
%%%%%%%                7.Minimum central pressure (mb)
%%%%%%%                8.Landfall
%%%%%%%                9-12.R34 in different quadtants
zz=find(isnan(Typhoon_data(:,3)));
Typhoon_data(zz,:)=[];
Ty=Typhoon_data;
%%%%% exclude the invalid values
zz=find(isnan(Ty(:,5))|isnan(Ty(:,7)));
Ty(zz,:)=[];
clearvars Typhoon_data
%%%%%% extract different TC tracjeotries
rank(1,1)=Ty(1,1);
rank(1,2)=1;
n=1;
for i=1:size(Ty,1)-1
    if Ty(i+1,1)~=Ty(i,1)
        rank(n,3)=i;
        rank(n+1,1)=Ty(i+1,1);
        rank(n+1,2)=i+1;
        n=n+1;
    end
end
rank(end,3)=size(Ty,1);
rank(:,4)=rank(:,3)-rank(:,2)+1;
%%%%%% calculate the move velocity of TC
for i=1:size(rank,1)
    for j=1:rank(i,4)
        if j==1
            Vt_u=(Ty(rank(i,2)-1+2,3)-Ty(rank(i,2),3))*111176*cos(Ty(rank(i,2)-1+2,4)-Ty(rank(i,2),4)/2*pi/180)/...
                ((Ty(rank(i,2)-1+2,2)-Ty(rank(i,2),2))*3600*24);%%%% 
            Vt_v=(Ty(rank(i,2)-1+2,4)-Ty(rank(i,2),4))*111176/((Ty(rank(i,2)-1+2,2)-Ty(rank(i,2),2))*3600*24);
        else
            Vt_u=(Ty(rank(i,2)-1+j,3)-Ty(rank(i,2)-2+j,3))*111176*cos(Ty(rank(i,2)-1+j,4)-Ty(rank(i,2)-2+j,4)/2*pi/180)/...
                ((Ty(rank(i,2)-1+j,2)-Ty(rank(i,2)-2+j,2))*3600*24);%%%%
            Vt_v=(Ty(rank(i,2)-1+j,4)-Ty(rank(i,2)-2+j,4))*111176/((Ty(rank(i,2)-1+j,2)-Ty(rank(i,2)-2+j,2))*3600*24);
        end
        Ty(rank(i,2)-1+j,13:14)=[Vt_u,Vt_v];
    end
end

%%%%%%%%%%%interpolation the TC data to hourly intervel 
Ty_data=cell(n,1);
for i=1:n  
    Ty1=Ty(rank(i,2):rank(i,3),:);
    if size(Ty1,1)>1
    tt=[Ty1(1,2):1/24:Ty1(end,2)];
    time=datevec(Ty1(:,2));
    Tyd=nan(length(tt),14);
    for j=1:rank(i,4)
        kk=find(tt==Ty1(j,2));
        if isempty(kk)==0
            Tyd(kk,:)=Ty1(j,:);
        end
    end
    Tyd(:,1)=fillmissing(Tyd(:,1),'linear');   
    for kk=1:2
        Tyd(:,kk)=fillmissing(Tyd(:,kk),'linear');
    end
    for kk=3:size(Tyd,2)
        Tyd(:,kk)=interp1(Ty1(:,2),Ty1(:,kk),Tyd(:,2));
    end
    else 
        Tyd=Ty1;
    end
    Ty_data{i,1}=Tyd;
end
Ty_d=cell2mat(Ty_data);
zz=find(Ty_d(:,6)<0);
[Ty_d,I1]=sortrows(Ty_d,2);
time=datevec(Ty_d(:,2));
%%%%%%% extract different TC tracjeotries again
ran(1,1)=time(1,2);
ran(1,2)=1;
nn=1;
for i=1:size(time,1)-1
    if time(i,2)~=time(i+1,2)
        ran(nn,3)=i;
        ran(nn+1,1)=time(i+1,2);
        ran(nn+1,2)=i+1;
        nn=nn+1;
    end
end
ran(end,3)=size(time,1);
ran(:,4)=ran(:,3)-ran(:,2)+1;
%%%%%%%%% read ERA5 data
for i=1:size(ran,1)
    i
    Data=Ty_d(ran(i,2):ran(i,3),:);
    time=datevec(Data(:,2));
    year=time(1,1);month=time(1,2);
    if month<10
        Month=[num2str(0) num2str(month)];
    else
        Month=num2str(month);
    end
    extra='V:\wind\ERA5\' ; %%% file path
    filename=['ERA199301.nc'];
    lat=ncread([extra,filename],'latitude');
    lon=ncread([extra,filename],'longitude');
    filename=['ERA' num2str(year) Month '.nc'];
    U_day=ncread([extra,filename],'u10');
    V_day=ncread([extra,filename],'v10');
    ZZ=nan(ran(i,4),1);
    for j=1:ran(i,4)      
        day=time(j,3);
        hour=time(j,4);
        %%%%%  Extraction of map data for plus or minus 6 degrees of latitude 
        %%%%%  and longitude of the typhoon center
        lon1=Data(j,3)-6;lon2=Data(j,3)+6;
        lat1=Data(j,4)-6;lat2=Data(j,4)+6;
        %%%%%% Combine maps in case of typhoons spanning 0 longitude
        if lon1>=0&lon2<=360
            lon_in=find((lon>=lon1)&(lon<=lon2));
            lat_in=find((lat>=lat1)&(lat<=lat2));
            lon3=double(lon(lon_in));
            lat3=double(lat(lat_in));
            u_day=U_day(lon_in,lat_in,24*(day-1)+hour+1);
            v_day=V_day(lon_in,lat_in,24*(day-1)+hour+1);
        elseif lon1<0
            lon_in1=find((lon>=lon1)&(lon<=lon2));
            lat_in=find((lat>=lat1)&(lat<=lat2));
            lon3a=double(lon(lon_in1));
            lat3=double(lat(lat_in));
            u_day1=U_day(lon_in1,lat_in,24*(day-1)+hour+1);
            v_day1=V_day(lon_in1,lat_in,24*(day-1)+hour+1);
            lon1b=lon1+360;
            lon_in2=find((lon>=lon1b)&(lon<=360));
            lon3b=double(lon(lon_in2));
            u_day2=U_day(lon_in2,lat_in,24*(day-1)+hour+1);
            v_day2=V_day(lon_in2,lat_in,24*(day-1)+hour+1);
            u_day=[u_day2;u_day1]; v_day=[u_day2;v_day1];
            lon3=[lon3b-360;lon3a];
        elseif lon2>360
            lon_in1=find((lon>=lon1)&(lon<=lon2));
            lat_in=find((lat>=lat1)&(lat<=lat2));
            lon3a=double(lon(lon_in1));
            lat3=double(lat(lat_in));
            u_day1=U_day(lon_in1,lat_in,24*(day-1)+hour+1);
            v_day1=V_day(lon_in1,lat_in,24*(day-1)+hour+1);
            lon2b=lon2-360;
            lon_in2=find((lon>=0)&(lon<=lon2b));
            lon3b=double(lon(lon_in2));
            u_day2=U_day(lon_in2,lat_in,24*(day-1)+hour+1);
            v_day2=V_day(lon_in2,lat_in,24*(day-1)+hour+1);
            u_day=[u_day1;u_day2]; v_day=[u_day1;v_day2];
            lon3=[lon3a;lon3b+360];
        end
        %%%%%%The left and right flip matrices, since the latitude is given by 90->-90
        lat3=flipud(lat3);u_day=fliplr(u_day);v_day=fliplr(v_day);
        vel=sqrt(u_day.^2+v_day.^2);
        vor=ra_vorticity(lat3,lon3,u_day,v_day);
        %%%%%%For the southern hemisphere, reverse it.
        zz=find(lat3<0);
        vor(:,zz)=-vor(:,zz);
        str=ra_strain_rate(lat,lon,u_day,v_day);
        %%%%% Refine the position of TC in ERA5
        Lon=Data(j,3); Lat=Data(j,4);
        Lon1=Lon-1;Lon2=Lon+1;
        Lat1=Lat-1;Lat2=Lat+1;
        Lon_in=find((lon3>=Lon1)&(lon3<=Lon2));
        Lat_in=find((lat3>=Lat1)&(lat3<=Lat2));
        Lon_t=lon3(Lon_in);Lat_t=lat3(Lat_in);
        Str=str(Lon_in,Lat_in);
        Vel=vel(Lon_in,Lat_in);
        Vor=vor(Lon_in,Lat_in);
        [a,b]=find(Vel==min(Vel(:))&Vor>=nanmean(Vor(:)));
        V_min=Vel(a,b);
        %%%%Ensure the accurate position whose wind speed is minimized and
        %%%% relative vorticity is more than the regional mean
        if length(a)==1
            Lon_c=Lon_t(a);Lat_c=Lat_t(b);
        elseif length(a)>1
            [a1,b1]=min(abs(a-size(Lon_t,1)/2)+abs(b-size(Lat_t,1)/2));
            Lon_c=Lon_t(a(b1));Lat_c=Lat_t(b(b1));
        elseif isempty(a)==1
        %%%%% use the original center position if you don't find the right conditions    
            Lon_c=Lon;Lat_c=Lat; 
        end
        [LAT,LON]=meshgrid(lat3,lon3);
        %%%%%%%%=====================================================================
        %%%%% Construct the wind field
        %%%%%%%%=====================================================================
        Pn=1.013*10^5; %%% standard atmosphere pressure，unit:Pa
        Pc=Data(j,7)*100; %%%% central atmosphere pressure ,unit:mb-->Pa
        Vmax=Data(j,5); %%%% max wind speed,unit:m/s
        %%%%%%% Radius of max wind speed，unit:nmile-->km，if the RMW is invalid, use the
        %%%%%%% formula to attain
        if isnan(Data(j,6))==0
            Rmax=Data(j,6)*1.852;
        else
            Rmax=51.6*exp(-0.0223*Vmax+0.0281*abs(Data(j,4)));
        end
        B=1.0036+0.0173*Vmax-0.0313*log(Rmax)+0.0087*abs(Data(j,4)); %%%% parameter B
        A=Rmax^B; %%%% parameter A
        rho=1.3; %%%%% air density:kg/m3
        w=2*pi/86400;
        f=sin(Data(j,4)*pi/180)*2*w;
        Vt_u=Data(j,13);Vt_v=Data(j,14);
        %%%%%%% ==================
        [~,z1]=min(abs(lon3-Lon_c));[~,z2]=min(abs(lat3-Lat_c));
           lon_c=lon3(z1);
        lat_c=lat3(z2);
        %%%%combine the models============
        dlon=LON-lon_c;
        dlat=LAT-lat_c;
        sita=nan(size(dlon));
        %%%%%%Willoughby adjustment method
        u_x=nan(size(dlon));
        v_y=nan(size(dlon));
        V_u_W=nan(size(dlon));
        V_v_W=nan(size(dlon));
        for ii=1:length(lon3)
            for jj=1:length(lat3)
                sita(ii,jj)=atan(dlat(ii,jj)/dlon(ii,jj));
                if dlon(ii,jj)<0
                    sita(ii,jj)=sita(ii,jj)+pi;
                elseif dlon(ii,jj)>0&dlat(ii,jj)<0
                    sita(ii,jj)=sita(ii,jj)+2*pi;
                end
                r=sqrt(dlon(ii,jj).^2+dlat(ii,jj).^2)*111;%%%%%unit: km
                Vg=Vmax.*((Rmax./r).^B.*exp(1-(Rmax./r).^B)).^0.5;
                if lat_c>=0                
                u_x(ii,jj)=-Vg*sin(sita(ii,jj));
                v_y(ii,jj)=Vg*cos(sita(ii,jj));
                else  
                %%%% The southern hemisphere rotates in the opposite direction
                u_x(ii,jj)=Vg*sin(sita(ii,jj));
                v_y(ii,jj)=-Vg*cos(sita(ii,jj));
                end
                c=r/(4*Rmax);%%%% parameter:n=4
                afa=c.^4/(1+c.^4);%%%% parameter:α
                if Vmax>33
                    V_u_W(ii,jj)=(1-afa).*(u_x(ii,jj)+Vt_u)+afa*u_day(ii,jj);
                    V_v_W(ii,jj)=(1-afa).*(v_y(ii,jj)+Vt_v)+afa*v_day(ii,jj);
                else
                    V_u_W(ii,jj)=(1-afa).*(u_x(ii,jj))+afa*u_day(ii,jj);
                    V_v_W(ii,jj)=(1-afa).*(v_y(ii,jj))+afa*v_day(ii,jj);
                end
            end
        end
        VW=sqrt(u_x.^2+v_y.^2);
        %%%%%%%Li adjustment method
        V_uL=nan(size(dlon));
        V_vL=nan(size(dlon));
        vel_e=vel(z1-8:z1+8,z2-8:z2+8);
        vel_e=sort(vel_e(:),'descend');
        vel_m=mean(vel_e(1:floor(length(vel_e)*0.05)));
        ratio=Vmax/vel_m;
        for ii=1:length(lon3)
            for jj=1:length(lat3)
                r=sqrt(dlon(ii,jj).^2+dlat(ii,jj).^2)*111;
                if r<Rmax
                    V_uL(ii,jj)=(r/Rmax*ratio+(Rmax-r)/Rmax)*u_day(ii,jj);
                    V_vL(ii,jj)=(r/Rmax*ratio+(Rmax-r)/Rmax)*v_day(ii,jj);
                elseif r>=Rmax&r<4*Rmax
                    V_uL(ii,jj)=((r-Rmax)/(3*Rmax)+(4*Rmax-r)/(3*Rmax)*ratio)*u_day(ii,jj);
                    V_vL(ii,jj)=((r-Rmax)/(3*Rmax)+(4*Rmax-r)/(3*Rmax)*ratio)*v_day(ii,jj);
                else
                    V_uL(ii,jj)=u_day(ii,jj);
                    V_vL(ii,jj)=v_day(ii,jj);
                end
            end
        end
%%%%%%% combine two models
        omega=(Data(j,8)-0)/(100-0);
        if omega<0
            omega=0;
        elseif omega>1
            omega=1;
        end
        V_u=omega*V_u_W+(1-omega)*V_uL;
        V_v=omega*V_v_W+(1-omega)*V_vL;

        V_u(z1,z2)=u_day(z1,z2);V_v(z1,z2)=v_day(z1,z2);
        V_com=sqrt(V_u.^2+V_v.^2);  
        %%%%%%%%%%%=====insect the corrected wind field to initial ERA5=========
        V_u=fliplr(double(V_u));V_v=fliplr(double(V_v));
         if lon1>=0&lon2<=360
            lon_in=find((lon>=lon1)&(lon<=lon2));
            lat_in=find((lat>=lat1)&(lat<=lat2));
            U_day(lon_in,lat_in,24*(day-1)+hour+1)=V_u;
            V_day(lon_in,lat_in,24*(day-1)+hour+1)=V_v;
        elseif lon1<0
            lon_in1=find((lon>=lon1)&(lon<=lon2));
            lat_in=find((lat>=lat1)&(lat<=lat2));
            U_day(lon_in1,lat_in,24*(day-1)+hour+1)=V_u(end-length(lon_in1)+1:end,:);
            V_day(lon_in1,lat_in,24*(day-1)+hour+1)=V_v(end-length(lon_in1)+1:end,:);
            lon1b=lon1+360;
            lon_in2=find((lon>=lon1b)&(lon<=360));
            U_day(lon_in2,lat_in,24*(day-1)+hour+1)=V_u(1:length(lon_in2),:);
            V_day(lon_in2,lat_in,24*(day-1)+hour+1)=V_v(1:length(lon_in2),:);
        elseif lon2>360
            lon_in1=find((lon>=lon1)&(lon<=lon2));
            lat_in=find((lat>=lat1)&(lat<=lat2));
            U_day(lon_in1,lat_in,24*(day-1)+hour+1)=V_u(1:length(lon_in1),:);
            V_day(lon_in1,lat_in,24*(day-1)+hour+1)=V_v(1:length(lon_in1),:);
            lon2b=lon2-360;
            lon_in2=find((lon>=0)&(lon<=lon2b));
            U_day(lon_in2,lat_in,24*(day-1)+hour+1)=V_u(end-length(lon_in2)+1:end,:);
            V_day(lon_in2,lat_in,24*(day-1)+hour+1)=V_v(end-length(lon_in2)+1:end,:);
         end
    end
    %% %%%%% create new nc.file
        MA=['G:\ERA5_RE_center\ERA'  num2str(year) Month '.nc'];     
        nccreate(MA,'lon','Dimensions',{'dimdx',size(U_day,1)});
        nccreate(MA,'lat','Dimensions',{'dimdy',size(U_day,2)});
        nccreate(MA,'u10','Dimensions',{'dimdx',size(U_day,1),'dimdy',size(U_day,2),'dimdt',size(U_day,3)},'Datatype','double');
        nccreate(MA,'v10','Dimensions',{'dimdx',size(U_day,1),'dimdy',size(U_day,2),'dimdt',size(U_day,3)},'Datatype','double');
        ncwrite(MA,'lon',lon);ncwriteatt(MA,'lon','units','degree_east');
        ncwrite(MA,'lat',lat);ncwriteatt(MA,'lat','units','degree_north');
        ncwrite(MA,'u10',U_day);ncwriteatt(MA,'u10','units','m/s');
        ncwrite(MA,'v10',V_day);ncwriteatt(MA,'v10','units','m/s');
        close
end