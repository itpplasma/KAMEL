%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = plot_dump(name, show_vertices)

% function parses the dump file, plots the partition of the region and also zeros found by ZeroSolver
% name - name of the dump file
% show_vertices - if it is positive then the edges vertices are additionally shown (may be slow for large files)

file = fopen(name, 'r');

skip_to_line('Re{zero}', file); % seek the file position to the first line containing substring 'Re{zero}'

zeros = get_list_of_vertices(file); % get list of zeros from the file

R = {}; % array of rectangles

n = 0;
while(true) % loop over rectangles

    if (show_vertices) % case when the edges vertices are needed

        line = skip_to_line('The rectangle winding number', file);

        if(~ischar(line)) break; end % check for end of the file

        n = n + 1; % increase the rectangles counter

        skip_to_line('Re{vertex}', file);
        R{n}.edge1 = get_list_of_vertices(file);

        skip_to_line('Re{vertex}', file);
        R{n}.edge2 = get_list_of_vertices(file);

        skip_to_line('Re{vertex}', file);
        R{n}.edge3 = get_list_of_vertices(file);

        skip_to_line('Re{vertex}', file);
        R{n}.edge4 = get_list_of_vertices(file);

        disp(['rectangle # ' num2str(n) ' processed successfully...']);

    else

        skip_to_line('The rectangle vertices:', file);

        line = skip_to_line('Re{vertex}', file);

        if(~ischar(line)) break; end % check for end of the file

        n = n + 1; % increase the rectangles counter

        corner = get_list_of_vertices(file); % the corner vertices of the rectangle

        R{n}.edge1.ind = [corner.ind(1) corner.ind(2)];
        R{n}.edge1.rez = [corner.rez(1) corner.rez(2)];
        R{n}.edge1.imz = [corner.imz(1) corner.imz(2)];
        R{n}.edge1.ref = [corner.ref(1) corner.ref(2)];
        R{n}.edge1.imf = [corner.imf(1) corner.imf(2)];
        R{n}.edge1.arg = [corner.arg(1) corner.arg(2)];

        R{n}.edge2.ind = [corner.ind(2) corner.ind(3)];
        R{n}.edge2.rez = [corner.rez(2) corner.rez(3)];
        R{n}.edge2.imz = [corner.imz(2) corner.imz(3)];
        R{n}.edge2.ref = [corner.ref(2) corner.ref(3)];
        R{n}.edge2.imf = [corner.imf(2) corner.imf(3)];
        R{n}.edge2.arg = [corner.arg(2) corner.arg(3)];

        R{n}.edge3.ind = [corner.ind(4) corner.ind(3)];
        R{n}.edge3.rez = [corner.rez(4) corner.rez(3)];
        R{n}.edge3.imz = [corner.imz(4) corner.imz(3)];
        R{n}.edge3.ref = [corner.ref(4) corner.ref(3)];
        R{n}.edge3.imf = [corner.imf(4) corner.imf(3)];
        R{n}.edge3.arg = [corner.arg(4) corner.arg(3)];

        R{n}.edge4.ind = [corner.ind(1) corner.ind(4)];
        R{n}.edge4.rez = [corner.rez(1) corner.rez(4)];
        R{n}.edge4.imz = [corner.imz(1) corner.imz(4)];
        R{n}.edge4.ref = [corner.ref(1) corner.ref(4)];
        R{n}.edge4.imf = [corner.imf(1) corner.imf(4)];
        R{n}.edge4.arg = [corner.arg(1) corner.arg(4)];

        disp(['rectangle # ' num2str(n) ' processed successfully...']);

    end

end

figure('Name', 'The dump of ZeroSolver data'); box on; hold on;

plot(zeros.rez, zeros.imz, 'ro');

cols = ['r', 'g', 'b', 'k', 'c', 'm', 'y'];

for k = 1:n % loop over rectangles

    col = cols(mod(k-1, 7)+1);

    plot(R{k}.edge1.rez, R{k}.edge1.imz, [col '.'], R{k}.edge1.rez, R{k}.edge1.imz, [col '-']);
    plot(R{k}.edge2.rez, R{k}.edge2.imz, [col '.'], R{k}.edge2.rez, R{k}.edge2.imz, [col '-']);
    plot(R{k}.edge3.rez, R{k}.edge3.imz, [col '.'], R{k}.edge3.rez, R{k}.edge3.imz, [col '-']);
    plot(R{k}.edge4.rez, R{k}.edge4.imz, [col '.'], R{k}.edge4.rez, R{k}.edge4.imz, [col '-']);
end

axis tight;
xlabel('Re\{z\}');
ylabel('Im\{z\}');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = skip_to_line(str, file)

% function reads the 'file' from the current position line by line until substring 'str' is found
% str - string to be found
% file - the file to be read
% return - the line where substring 'str' is found

line = '';

while(ischar(line))

    line = fgetl(file);

    if (~isempty(strfind(line, str))) break; end

end

w = line;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function list = get_list_of_vertices(file)

% function reads the 'file' from the current position line by line until an empty line is found
% file - the file to be read
% return - a structure that holds list of vertices where:
% ind - vertex index
% rez - real part of z
% imz - imag part of z
% ref - real part of f(z)
% imf - imag part of f(z)
% arg - argument of f(z)

list.ind = 0;
list.rez = 0;
list.imz = 0;
list.ref = 0;
list.imf = 0;
list.arg = 0;

n = 1;

while(true)

    line = fgetl(file);

    if (~length(line)) break; end % check for empty line

    remain = line;

    % parse line to extract ind, z, f(z) as written by ZeroSolver:
    [str, remain] = strtok(remain); list.ind(n) = str2num(str);
    [str, remain] = strtok(remain); list.rez(n) = str2num(str);
    [str, remain] = strtok(remain); list.imz(n) = str2num(str);
    [str, remain] = strtok(remain); list.ref(n) = str2num(str);
    [str, remain] = strtok(remain); list.imf(n) = str2num(str);
    [str, remain] = strtok(remain); list.arg(n) = str2num(str);

    n = n+1;

end

list.ind = list.ind';
list.rez = list.rez';
list.imz = list.imz';
list.ref = list.ref';
list.imf = list.imf';
list.arg = list.arg';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
