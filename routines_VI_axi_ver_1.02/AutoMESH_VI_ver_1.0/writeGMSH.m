function [ fileID ] = writeGMSH(filename,GEOMETRY,SPACING)

if abs(GEOMETRY(1,1)-GEOMETRY(end,1)) < 1E-10 && ...
     abs(GEOMETRY(1,2)-GEOMETRY(end,2)) < 1E-10
    GEOMETRY=GEOMETRY(1:end-1,:);
end

fileID = fopen([filename '.geo'],'w');

shape{1}=GEOMETRY;
element=shape{1};
spacing=zeros(size(element,1),1);
spacing(:)=SPACING;

%% POINTS
n_point=1;
for ii = 1:size(element,1)
    fprintf(fileID,'%s  %d %s %8f %s %8f %s %8f %s \n',...
        char('Point('),...
        n_point,...
        char(')={'),...
        element(ii,1),...
        ',',...
        element(ii,2),...
        ',0,',...
        spacing(ii),...
        '};');
    n_point=n_point+1;
end

fprintf(fileID,'%6s \n',' ');

%% LINES
n_line=1;
for ii = 1:size(element,1)
    if ii<size(element,1)
        fprintf(fileID,'%s  %d %s %d %s %d %s \n',...
            char('Line('),...
            n_line,...
            char(')={'),...
            n_line,...
            ',',...
            n_line+1,...
            '};');
    elseif ii==size(element,1)
        fprintf(fileID,'%s  %d %s %d %s %d %s \n',...
            char('Line('),...
            n_line,...
            char(')={'),...
            n_line,...
            ',',...
            n_line-size(element,1)+1,...
            '};');
    end
    n_line=n_line+1;
end

fprintf(fileID,'%1s \n',' ');

%% LINE LOOP
n_loop=1;
fprintf(fileID,'%s \n','Line Loop(1)={');
for ii = 1:size(element,1)
    if ii<size(element,1)
        fprintf(fileID,'%d %s\n',...
            n_line-size(element,1)+ii-1,...
            ',');
    elseif ii==size(element,1)
        fprintf(fileID,'%d %s\n',...
            n_line-size(element,1)+ii-1, ...
            '};');
    end
end

fprintf(fileID,'%1s \n',' ');

%% PLANE SURFACE
fprintf(fileID,'%s \n','Plane Surface(1)=(1);');


fclose(fileID);

end

