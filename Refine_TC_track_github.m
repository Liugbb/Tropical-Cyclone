clc
clear
load('V:\Typhoon\IBTrACS\Typhoon_data_new')
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
%%%%%%%%===================================================================
zz=find(isnan(Typhoon_data(:,3)));
Typhoon_data(zz,:)=[];
Ty=Typhoon_data;
zz=find(isnan(Ty(:,5))|isnan(Ty(:,7)));%%%%%
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
                ((Ty(rank(i,2)-1+2,2)-Ty(rank(i,2),2))*3600*24);%%%% 考虑纬度收缩
            Vt_v=(Ty(rank(i,2)-1+2,4)-Ty(rank(i,2),4))*111176/((Ty(rank(i,2)-1+2,2)-Ty(rank(i,2),2))*3600*24);
        else
            Vt_u=(Ty(rank(i,2)-1+j,3)-Ty(rank(i,2)-2+j,3))*111176*cos(Ty(rank(i,2)-1+j,4)-Ty(rank(i,2)-2+j,4)/2*pi/180)/...
                ((Ty(rank(i,2)-1+j,2)-Ty(rank(i,2)-2+j,2))*3600*24);%%%% 考虑纬度收缩
            Vt_v=(Ty(rank(i,2)-1+j,4)-Ty(rank(i,2)-2+j,4))*111176/((Ty(rank(i,2)-1+j,2)-Ty(rank(i,2)-2+j,2))*3600*24);
        end
        Ty(rank(i,2)-1+j,13:14)=[Vt_u,Vt_v];
    end
end

% % % [Ty_d,I1]=sortrows(Ty,2);
% % % time=datevec(Ty_d(:,2));
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
     %%%%%% For close to land for nan values, the principle of proximity interpolation is used
    Tyd(:,8)=fillmissing(Tyd(:,8),'nearest');
    Ty_data{i,1}=Tyd;
end
Ty_d=cell2mat(Ty_data);

zz=find(Ty_d(:,6)<0);
[Ty_d,I1]=sortrows(Ty_d,2);
time=datevec(Ty_d(:,2));
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
 TCs_BT=cell(size(ran,1),1);
for i=1:314
    i
    Data=Ty_d(ran(i,2):ran(i,3),:);
    time=datevec(Data(:,2));
    year=time(1,1);month=time(1,2);
    if month<10
        Month=[num2str(0) num2str(month)];
    else
        Month=num2str(month);
    end
    extra='V:\wind\ERA5\' ;
    filename=['ERA199301.nc'];
    lat=ncread([extra,filename],'latitude');
    lon=ncread([extra,filename],'longitude');
    extra1='F:\ERA5_RE_center_new_double\' ;
    filename1=['ERA' num2str(year) Month '.nc'];
    U_day=ncread([extra1,filename1],'u10');
% % %     UUU=U_day(:,:,600);
    V_day=ncread([extra1,filename1],'v10');
    ZZ=nan(ran(i,4),1);
    TCs=nan(size(Data,1,10));
    for j=1:ran(i,4)  
        day=time(j,3);
        hour=time(j,4);
        %%%%%  Extraction of map data for plus or minus 8 degrees of latitude 
        %%%%%  and longitude of the typhoon center
        lon1=Data(j,3)-8;lon2=Data(j,3)+8;
        lat1=Data(j,4)-8;lat2=Data(j,4)+8;
        %%%%%%======遇到跨0经度的台风情况=================
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
        %%%% For the southern hemisphere, reverse it.
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
            Lon_c=Lon;Lat_c=Lat;   %%%%没有找到合适条件的话就用原来的涡心位置
        end
        [LAT,LON]=meshgrid(lat3,lon3);
        %%%%%%%%=====================================================================
        %%%%% 对数据进行修正处理
        %%%%%%%%=====================================================================
        Pn=1.013*10^5; %%% standard atmosphere pressure，unit:Pa
        Pc=Data(j,7)*100; %%%%central atmosphere pressure,unit:mb-->Pa
        Vmax=Data(j,5); %%%% max wind speed,unit:m/s
        %%%%% use the original center position if you don't find the right conditions    
         if isnan(Data(j,6))==0
            Rmax=Data(j,6)*1.852;
        else
            Rmax=51.6*exp(-0.0223*Vmax+0.0281*abs(Data(j,4)));
        end
        B=1.0036+0.0173*Vmax-0.0313*log(Rmax)+0.0087*abs(Data(j,4)); %%%% parameter B
        A=Rmax^B; %%%% parameter A
        rho=1.3; %%%%% air density, unit:kg/m3
        w=2*pi/86400;
        f=sin(Data(j,4)*pi/180)*2*w;%%%%% 科氏参数
        Vt_u=Data(j,13);Vt_v=Data(j,14);
        [~,z1]=min(abs(lon3-Lon_c));[~,z2]=min(abs(lat3-Lat_c));
        lon_c=lon3(z1);
        lat_c=lat3(z2);
        %%%%%%% calculate max wind spped in ERA5=============================
        vel_e=vel(z1-8:z1+8,z2-8:z2+8);
        vel_e=sort(vel_e(:),'descend');
        %%%% Find the maximum threshold for comparison with IBTrACS.
        vel_m=mean(vel_e(1:floor(length(vel_e)*0.04)));
        %%%%模型融合============
        dlon=LON-lon_c;
        dlat=LAT-lat_c;
        ddis=sqrt(dlon.^2+dlat.^2);
        

        %%%%取一个纬向剖面
        ddis=sqrt(dlon.^2+dlat.^2);
        dis=0:0.25:10;
        vel_r=nan(length(dis),2);
        for m=1:length(dis)
            
            [a,b]=find(ddis>=dis(m)-0.25/2&ddis<=dis(m)+0.25/2);
            vv=[];
            for n=1:length(a)
                vv(n)=vel(a(n),b(n));
            end
            vel_r(m,1)=dis(m);
            vel_r(m,2)=nanmean(vv);
        end

        [Vmax_RE,RMW_in]=max(vel_r(:,2));  %%%% MWS
        RMW_RE=vel_r(RMW_in,1)*111;    %%%% RMW
        
        %%%%提取七级风圈半径
        %%%%%% NE quadrant
        dis=0:0.25:10;
        vel_NE=nan(length(dis),2);
        for m=1:length(dis)
              
            [a,b]=find(ddis>=dis(m)-0.25/2&ddis<=dis(m)+0.25/2&...
                dlon>=0&dlat>=0);
            vv=[];
            for n=1:length(a)
                vv(n)=vel(a(n),b(n));
            end
            vel_NE(m,1)=dis(m);
            vel_NE(m,2)=nanmean(vv);
        end
        vel_NE1=vel_NE(:,2)-17; 
        vel_NE2=(vel_NE1(1:end-1).*vel_NE1(2:end));
        zz1=find(vel_NE2<0);
        zz1(zz1<RMW_in)=[];
        if isempty(zz1)==0
        R34_NE=(vel_NE(zz1(1),1)+vel_NE(zz1(1)+1,1))/2*111;
        else
        R34_NE=nan;
        end
        %%%%%% SE quadrant
        dis=0:0.25:10;
        vel_SE=nan(length(dis),2);
        for m=1:length(dis)
            
            [a,b]=find(ddis>=dis(m)-0.25/2&ddis<=dis(m)+0.25/2&...
                dlon>=0&dlat<=0);
            vv=[];
            for n=1:length(a)
                vv(n)=vel(a(n),b(n));
            end
            vel_SE(m,1)=dis(m);
            vel_SE(m,2)=nanmean(vv);
        end
        vel_SE1=vel_SE(:,2)-17; 
        vel_SE2=(vel_SE1(1:end-1).*vel_SE1(2:end));
        zz2=find(vel_SE2<0);
        zz2(zz2<RMW_in)=[];
       if isempty(zz2)==0
        R34_SE=(vel_SE(zz2(1),1)+vel_SE(zz2(1)+1,1))/2*111;
        else
        R34_SE=nan;
        end
        %%%%%% SW quadrant
        dis=0:0.25:10;
        vel_SW=nan(length(dis),2);
        for m=1:length(dis)
            
            [a,b]=find(ddis>=dis(m)-0.25/2&ddis<=dis(m)+0.25/2&...
                dlon<=0&dlat<=0);
            vv=[];
            for n=1:length(a)
                vv(n)=vel(a(n),b(n));
            end
            vel_SW(m,1)=dis(m);
            vel_SW(m,2)=nanmean(vv);
        end
        vel_SW1=vel_SW(:,2)-17; 
        vel_SW2=(vel_SW1(1:end-1).*vel_SW1(2:end));
        zz3=find(vel_SW2<0);
        zz3(zz3<RMW_in)=[];
        if isempty(zz3)==0
        R34_SW=(vel_SW(zz3(1),1)+vel_SW(zz3(1)+1,1))/2*111;
        else
        R34_SW=nan;
        end
        %%%%%% NW quadrant
        dis=0:0.25:10;
        vel_NW=nan(length(dis),2);
        for m=1:length(dis)
            
            [a,b]=find(ddis>=dis(m)-0.25/2&ddis<=dis(m)+0.25/2&...
                dlon<=0&dlat>=0);
            vv=[];
            for n=1:length(a)
                vv(n)=vel(a(n),b(n));
            end
            vel_NW(m,1)=dis(m);
            vel_NW(m,2)=nanmean(vv);
        end
        vel_NW1=vel_NW(:,2)-17; 
        vel_NW2=(vel_NW1(1:end-1).*vel_NW1(2:end));
        zz4=find(vel_NW2<0);
        zz4(zz4<RMW_in)=[];
        if isempty(zz4)==0
        R34_NW=(vel_NW(zz4(1),1)+vel_NW(zz4(1)+1,1))/2*111;
        else
        R34_NW=nan;
        end
        %%%%%%%====================================================
        TCs(j,1:2)=Data(j,1:2);
        TCs(j,3)=Lon_c;
        TCs(j,4)=Lat_c;
        TCs(j,5)=Vmax_RE; %%%Vmax
        TCs(j,6)=RMW_RE; %%%Rmax
        TCs(j,7:10)=[R34_NE,R34_SE,R34_SW,R34_NW];
        TCs(j,11)=vel_m; %%%%V_maxRE
    end
    TCs_BT{i,1}=TCs;
end
% % % save TCs_BT TCs_BT  
TCs_BTD=cell2mat(TCs_BT);
TCs_BTD=sortrows(TCs_BTD,2);
TCs_BTD=sortrows(TCs_BTD,1);
save TCs_BTD_new TCs_BTD

