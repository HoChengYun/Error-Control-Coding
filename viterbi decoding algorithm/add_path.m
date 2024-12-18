function [path, distance] = add_path(path, distance, rec_signal)
    path_num = 0; level = 0; new_path = []; new_distance = [];
    temp_distance = 0; temp_output = logical([0,0]);
    rec_signal = reshape(rec_signal, [2,length(rec_signal)./2]);
    [path_num, level] = size(path);
    for i = 1:path_num
        switch path(i,end)
            case 0
                %s0->s0
                new_path = [new_path; path(i,:), 0]; %update new path
                temp_output = logical([0 0]);
                temp_distance = double(xor(temp_output(1), rec_signal(1, level)))...
                                    + double(xor(temp_output(2), rec_signal(2, level)));
                new_distance = [new_distance; distance(i) + temp_distance];
                %s0->s1
                temp_output = logical([1 1]);
                new_path = [new_path; path(i,:), 1];
                temp_distance = double(xor(temp_output(1), rec_signal(1, level)))...
                                    + double(xor(temp_output(2), rec_signal(2, level)));
                new_distance = [new_distance; distance(i) + temp_distance];
            case 1
                %s1->s2
                temp_output = logical([1 0]);
                new_path = [new_path; path(i,:), 2]; %update new path
                temp_distance = double(xor(temp_output(1), rec_signal(1, level)))...
                                    + double(xor(temp_output(2), rec_signal(2, level)));
                new_distance = [new_distance; distance(i) + temp_distance];
                %s1->s3
                temp_output = logical([0 1]);
                new_path = [new_path; path(i,:), 3];
                temp_distance = double(xor(temp_output(1), rec_signal(1, level)))...
                                    + double(xor(temp_output(2), rec_signal(2, level)));
                new_distance = [new_distance; distance(i) + temp_distance];
            case 2
                %s2->s0
                temp_output = logical([1 1]);
                new_path = [new_path; path(i,:), 0]; %update new path
                temp_distance = double(xor(temp_output(1), rec_signal(1, level)))...
                                    + double(xor(temp_output(2), rec_signal(2, level)));
                new_distance = [new_distance; distance(i) + temp_distance];
                %s2->s1
                temp_output = logical([0 0]);
                new_path = [new_path; path(i,:), 1];
                temp_distance = double(xor(temp_output(1), rec_signal(1, level)))...
                                    + double(xor(temp_output(2), rec_signal(2, level)));
                new_distance = [new_distance; distance(i) + temp_distance];
            case 3
                %s3->s2
                temp_output = logical([0 1]);
                new_path = [new_path; path(i,:), 2]; %update new path
                temp_distance = double(xor(temp_output(1), rec_signal(1, level)))...
                                    + double(xor(temp_output(2), rec_signal(2, level)));
                new_distance = [new_distance; distance(i) + temp_distance];
                %s3->s3
                temp_output = logical([1 0]);
                new_path = [new_path; path(i,:), 3];
                temp_distance = double(xor(temp_output(1), rec_signal(1, level)))...
                                    + double(xor(temp_output(2), rec_signal(2, level)));
                new_distance = [new_distance; distance(i) + temp_distance];
        end
    end
        %delete long distance path
        same_index = 0; min_index = 0; min_distance = 0;
        %clear distance & path
        distance = []; path = [];
        for check = 0:3
            same_index = find(new_path(:,end) == check);
            if length(same_index) == 1
                distance = [distance; new_distance(same_index)];
                path = [path; new_path(same_index,:)];
                continue
            end
            %min_index : index of same_index
            [min_distance, min_index] = min(new_distance([same_index]));
            distance = [distance; min_distance];
            path = [path; new_path(same_index(min_index),:)];
        end
end