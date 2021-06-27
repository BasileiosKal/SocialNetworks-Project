% Databases paths
DB_AR_PATH = "C:\Users\KALOUA\Desktop\Vasilhs\Databases\AR.db";
DB_EIU_PATH = "C:\Users\KALOUA\Desktop\Vasilhs\Databases\EIU.db";
DB_OVERLAPS_PATH = "C:\Users\KALOUA\Desktop\Vasilhs\Databases\Overlaps.db";
DB_USERS_PATH = "C:\Users\KALOUA\Desktop\Vasilhs\Databases\SubGraphNodes.db";
DB_RS_PATH = "C:\Users\KALOUA\Desktop\Vasilhs\Databases\RS_2v1.db";

% Global variables
global NUM_CATEGORIES
global NUM_USERS
NUM_USERS = 400;
NUM_CATEGORIES = 10;


% Connect to the databases
AR_conn = sqlite(DB_AR_PATH,'readonly');             % The average ratings
EIU_conn = sqlite(DB_EIU_PATH,'readonly');           % The reliable sets sizes
Overaps_conn = sqlite(DB_OVERLAPS_PATH,'readonly');  % The overlaps
Users_conn = sqlite(DB_USERS_PATH,'readonly');       % The Users with indexes
RS_conn = sqlite(DB_RS_PATH,'readonly');             % The reliable sets

disp("Conected to the databases")

%% Create matrices containing the constants of the LP/IP

% ---> Create a mapping between users ids and indexes on those matrices
Map_UsersToIndx = containers.Map;
% Data from the SubGraphNodes database
Users_data = fetch(Users_conn, "SELECT * FROM Nodes");

for row_i = 1:length(Users_data)
    user = Users_data(row_i, 1); 
    user = user{1};   % the data returned from fetch are of the form cell array
    index = Users_data(row_i, 2);
    index = index{1};
    Map_UsersToIndx(user) = index + 1;
end

disp("Done creating users to index maping")

% ---> Load the sizes of the reachable sets to a
%      NUM_CATEGORIES x NUM_USERS matrix
global Eiu_matrix;
Eiu_matrix = zeros(NUM_CATEGORIES, NUM_USERS);
for cat = 1:NUM_CATEGORIES
    % get reachable set sizes for that category
    Eiu_data = fetch(EIU_conn, sprintf("SELECT * FROM Lengths_%d",cat));
    
    for EIU_row = 1:length(Eiu_data)
        user = Eiu_data(EIU_row, 1);
        user = user{1};  % get the ellement from the cell array
        eiu = Eiu_data(EIU_row, 2);
        eiu = eiu{1};
        
        Eiu_matrix(cat, Map_UsersToIndx(user)) = eiu;
    end
end

disp("Done creating Eiu matrix")

% ---> Create the overlaps matrices
%      put the overlaps to NUM_CATEGORIES matricies 
%      of dimensions NUM_USERS x  NUM_USERS
global Overlaps_matrix
Overlaps_matrix = zeros(NUM_USERS, NUM_USERS, NUM_CATEGORIES);
for cat = 1:NUM_CATEGORIES
    Overlaps_data = fetch(Overaps_conn, sprintf("SELECT * FROM OV_%d", cat));
    
    for OV_row = 1:length(Overlaps_data)
        uid = Overlaps_data(OV_row, 1); uid = uid{1}; uid_indx = Map_UsersToIndx(uid);
        vid = Overlaps_data(OV_row, 2); vid = vid{1}; vid_indx = Map_UsersToIndx(vid);
        ov = Overlaps_data(OV_row, 3); ov = ov{1};
        Overlaps_matrix(uid_indx, vid_indx, cat) = ov;
    end
end

disp("Done creating overlaps matrix")

% ---> Create the Iu(i, v) (1 if user v  belongs to the
%      reliable set of u for category i, 0 otherwise)
global IUV_matrix
IUV_matrix = zeros(NUM_USERS, NUM_USERS, NUM_CATEGORIES);

for cat = 1:NUM_CATEGORIES
    RS_data = fetch(RS_conn, sprintf("SELECT * FROM RS_%d", cat));
    for RS_row = 1:length(RS_data)
        uid = RS_data(RS_row, 1); uid = uid{1}; uid_indx = Map_UsersToIndx(uid);
        vid = RS_data(RS_row, 2); vid = vid{1}; vid_indx = Map_UsersToIndx(vid);
        if uid_indx ~= vid_indx
            %disp("WTF")
            IUV_matrix(uid_indx, vid_indx, cat) = 1;
        end
    end
end

disp("Done creating the IUV matrix")


% ---> Create the AR matrix of sape NUM_CATEGORIES x NUM_USERS
global AR_matrix
AR_matrix = zeros(NUM_CATEGORIES, NUM_USERS);

for cat = 1:NUM_CATEGORIES
    AR_data = fetch(AR_conn, sprintf("SELECT * FROM AR_%d", cat));
    for AR_row = 1:length(AR_data)
        uid = AR_data(AR_row, 1); uid = uid{1}; uid_indx = Map_UsersToIndx(uid);
        ariu = AR_data(AR_row, 2); ariu = ariu{1};
        AR_matrix(cat, uid_indx) = ariu;
    end
end

disp("Done creating the AR matrix")

%% Create matrices corresponding to each constrain.

% create matrix for wich A*x(:) =< M will give the first constrain
Con1_matrix = zeros(NUM_USERS, NUM_CATEGORIES*NUM_USERS);

for user_i = 1:NUM_USERS
    i = user_i;
    for j = 1:NUM_CATEGORIES
        Con1_matrix(user_i, i) = 1;
        i = i + NUM_USERS;
    end
end


% create matrix for wich A*x(:) =< M' will give the second constrain
Con2_matrix = zeros(NUM_CATEGORIES, NUM_CATEGORIES*NUM_USERS);

for cat_i = 1:NUM_CATEGORIES
    i = 1 + (cat_i-1)*NUM_USERS;
    for j = 1:NUM_USERS
        Con2_matrix(cat_i, i) = 1;
        i = i + 1;
    end
end


% create matrix for wich A*x(:) =< L will give the third constrain
Con3_matrix = zeros(NUM_USERS, NUM_USERS*NUM_CATEGORIES);

for uid = 1:NUM_USERS
   col_i = 1;
   for vid = 1:NUM_USERS
       for cat = 1:NUM_CATEGORIES
          if uid ~= vid
              Con3_matrix(uid, col_i) = IUV_matrix(uid, vid, cat);
          else
              Con3_matrix(uid, col_i) = 1;
          end
          col_i = col_i + 1;
       end
   end
end

disp("Done creating the constrain matrices")
%% LP

cat_M= 60;
user_M = 10;
L = 10;
 
% Initial point
X_0 =  (1/10000)*ones(1, NUM_CATEGORIES*NUM_USERS);

% Create matrices A and b so that A*x(:) =< b will give all the constrains
A = [Con1_matrix; Con2_matrix; Con3_matrix];

b_con1 = user_M * ones(NUM_USERS, 1);
b_con2 = cat_M* ones(NUM_CATEGORIES, 1);
b_L = L * ones(NUM_USERS, 1);

b = [b_con1; b_con2; b_L];

% lb: lower bound for the solution, 
% ub: upper bound for the solution
lb = zeros(1, NUM_CATEGORIES * NUM_USERS);
ub = ones(1, NUM_CATEGORIES * NUM_USERS);
% Other parameters and options
Aeq = [];
beq = [];
nonlcon = [];
fun = @optimize;

fmincon_options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000, ...
    'Display', 'iter', 'StepTolerance', 1.0000e-5, 'FiniteDifferenceStepSize', sqrt(eps),...
    'ConstraintTolerance', 1.0000e-1, 'SpecifyObjectiveGradient',true, 'CheckGradients', false);

sopt_options = optimoptions(@surrogateopt,'MaxFunctionEvaluations',100000, ...
    'Display', 'iter',...
    'ConstraintTolerance', 1.0000e-1, 'InitialPoints', X_0);

% solving the LP
fprintf('Solving LP...')
% [surrogate_x,fval,exitflag,output] = surrogateopt(fun, lb, ub, [], A, b, [], [], sopt_options);
[x,fval,exitflag,output] = fmincon(fun, surrogate_x, A, b, Aeq, beq, lb, ub, [], fmincon_options);

fprintf('Done \n')

%% IP 
% ---> Get the integer solution from the LP solution

% vectors to keep track of the constrains
int_x2 = zeros(1, NUM_CATEGORIES*NUM_USERS);
con_1 = zeros(1, NUM_USERS);
con_2 = zeros(1, NUM_CATEGORIES);
con_3 = zeros(1, NUM_USERS);


for cat = 1:NUM_CATEGORIES
    % get the M largest elements from the solution for each category
    M_max_elements = maxk(x(1+(cat-1)*NUM_USERS: cat*NUM_USERS), cat_M);
    for indx = 1+(cat-1)*NUM_USERS: cat*NUM_USERS
        if mod(indx, NUM_USERS) > 0
            con_1_indx = mod(indx, NUM_USERS);
        else
            con_1_indx = NUM_USERS;
        end
        
        % Check the constrains and if there is not a put the user to the 
        % solution, else remove the user (by setting x(indx) = -1) and 
        % re-short the users.
        if x(indx) >= min(M_max_elements) && con_1(con_1_indx) <= user_M && con_3(con_1_indx) <= L
            int_x2(indx) = 1;
            con_1(con_1_indx) = con_1(con_1_indx) + 1;
            con_2(cat) = con_2(cat) + 1;
            con_3(con_1_indx) = con_3(con_1_indx) + Con3_matrix*int_x2(:);
        elseif con_1(con_1_indx) > user_M
            x(indx) = -1;
            M_max_elements = maxk(x(1+(cat-1)*NUM_USERS: cat*NUM_USERS), cat_M);
        end
    end
end

disp("Done creating int solution")
disp(fun(int_x2))


%% Ploting the opt function
% X = x;
fun = @optimize;
X = rand(1, NUM_CATEGORIES * NUM_USERS);
% [X_1, X_2] = meshgrid(0:0.05:1, 0:0.05:1);
dims_1 = 21;
dims_2 = 21;

X_1 = sort(rand(dims_1, 1));
X_2 = sort(rand(dims_2, 1));

Z = zeros(dims_1, dims_2);

ind_1 = 3;
ind_2 = 111;    % 5
for i = 1:dims_1
    disp(i)
    for j = 1:dims_2
        % Z(i, j) = -fun([X_1(i, j),X(2:410), X_2(i, j), X(412:end)]);
       %  X = (1/1000) * rand(1, NUM_CATEGORIES * NUM_USERS); 
        Z(i, j) = -fun([X(1:ind_1-1), X_1(i),X(ind_1+1:ind_2-1), X_2(j), X(ind_2+1:end)]);
    end
end
surf(X_1, X_2, Z)
%% close and clear connections to dbs
close(AR_conn)
close(EIU_conn)
close(Overaps_conn)
close(Users_conn)
close(RS_conn)
clear AR_conn
clear EIU_conn
clear EIU_conn
clear Users_conn
clear RS_conn

disp("Databases closed")


%% Functions
% ----------------------------------------------------------- %
% Function to maximize
function [y, Grad] = optimize(x)
global Eiu_matrix;
global AR_matrix;
global NUM_CATEGORIES;
global Overlaps_matrix;
global NUM_USERS

lambda_par = 0;
Grad = zeros(1, NUM_CATEGORIES*NUM_USERS);

% disp(size(x))
ov = 0;
starting_at = 1;
for cat = 1:NUM_CATEGORIES
    ov = ov + ...
        sum(Overlaps_matrix(1:NUM_USERS,1:NUM_USERS,cat).*...
        (x(starting_at:starting_at-1+NUM_USERS)' * x(starting_at:starting_at-1+NUM_USERS)), 'all');
    
    Grad(starting_at:starting_at-1+NUM_USERS) = (Overlaps_matrix(1:NUM_USERS,1:NUM_USERS,cat) *...
        (x(starting_at:starting_at-1+NUM_USERS)'))';
    
    starting_at = 1 + (cat-1)*NUM_USERS;
end

% y = 0.5*y - sum((Eiu_matrix.*x), 'all') - lambda_par * sum(AR_matrix.*x, 'all');
y_eiu = 0;
i = 1;

for cat = 1:NUM_CATEGORIES
    for uid = 1:NUM_USERS
        % if Eiu_matrix(cat, uid) ~= 0
        y_eiu = y_eiu + (Eiu_matrix(cat, uid) + (lambda_par*AR_matrix(cat, uid)))* x(i);
        % end
        Grad(i) = Grad(i) - Eiu_matrix(cat, uid) - (lambda_par*AR_matrix(cat, uid));
        i = i + 1;
    end
end
%y = -sum((Eiu_matrix.*x), 'all');
y = 0.5*ov - (y_eiu); % - (lambda_par * sum(AR_matrix.*x, 'all'));
end


% ----------------------------------------------------------- %
% Calculate the number of nodes reached by a solution for a category
function [nodesList, total_users] = nodes_reached(solution, category, RSdb_conn, user_to_indx)

RS_data2 = fetch(RSdb_conn, sprintf("SELECT DISTINCT(uid) FROM RS_%d", category));
% Number of users in the subgraph that have given a rating to the
% specified category
total_users = [];
for row_i = 1:length(RS_data2)
    uid_cell = RS_data2(1); user_id = uid_cell{1};
    if user_to_indx(user_id) <= length(solution)
        total_users = [total_users, sprintf("'%s'", user_id)];
        reached_nodes_data2 = fetch(RSdb_conn,...
                  sprintf("SELECT vid FROM RS_%d WHERE uid == '%s'", category, user_id));
        for vid_i = 1:length(reached_nodes_data2)
            vid_cell = reached_nodes_data2(vid_i); vid = vid_cell{1};
            str_vid = sprintf("'%s'", vid);
            if ~ismember(str_vid, total_users)
                total_users = [total_users, str_vid];
            end
        end
    end
end


RS_data = fetch(RSdb_conn, sprintf("SELECT DISTINCT(uid) FROM RS_%d", category));
nodesList = [];
for row_i = 1:length(RS_data)
    row_data = RS_data(row_i);

    uid = row_data{1};   % get the users id from the cell array returned from fetch
    uid_indx = user_to_indx(uid);
    if uid_indx <= length(solution)
        
        if solution(uid_indx) == 1
             % get the nodes reachable by uid
            reached_nodes_data = fetch(RSdb_conn,...
                  sprintf("SELECT vid FROM RS_%d WHERE uid == '%s'", category, uid));
           
              nodesList = [nodesList, sprintf("'%s'", uid)];
            
            for rn_i = 1:length(reached_nodes_data)
                
                reached_node_cell = reached_nodes_data(rn_i);
                reached_node = reached_node_cell{1};
                s_vid = sprintf("'%s'", reached_node);
                if ~ismember(s_vid, nodesList)
                    nodesList = [nodesList, s_vid];
                    % nodesList = nodesList + 1;
                end
            end
        end
    end
end
end

