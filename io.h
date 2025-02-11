#ifndef GNNS_IO_H
#define GNNS_IO_H

#include "util/matrix.h"
#include "general.h"
#include <string>

//namespace gnns
//{
    template <typename T>
    void save_to_file(const std::vector<std::vector<T>>& m, const std::string& file_name)
    {
        FILE *fp = fopen(file_name.c_str(), "wb");

        if(fp == NULL)
        {
            throw GnnsException(std::string("File ") + file_name + "is NOT exist!\n");
        }

        for(const auto& row:m)
        {
            //write the dim
            int cols=static_cast<int>(row.size());
            fwrite(&cols, 4, 1, fp);
            //write the row data
            fwrite(row.data(), sizeof(T)*cols, 1, fp);

        }

        fclose(fp);
    }

    void save_start_end_positions_to_file(const std::vector<std::vector<std::pair<int,int>>>& data,const std::string& file_name)
    {
        FILE *fp = fopen(file_name.c_str(), "wb");

        if (fp == nullptr)
        {
            throw std::runtime_error("Cannot open file: " + file_name);
        }

        int rows = static_cast<int>(data.size());
        fwrite(&rows, sizeof(int), 1, fp);

        for (int i = 0; i < rows; ++i)
        {
            int cols = static_cast<int>(data[i].size());
            fwrite(&cols, sizeof(int), 1, fp);

            for (int j = 0; j < cols; ++j)
            {
                int first = data[i][j].first;
                int second = data[i][j].second;
                fwrite(&first, sizeof(int), 1, fp);
                fwrite(&second, sizeof(int), 1, fp);
            }
        }

        fclose(fp);
    }

    void save_top_level_start_end_positions_to_file(const std::vector<std::vector<std::pair<int,int>>>& data,const std::string& file_name)
    {
        FILE *fp = fopen(file_name.c_str(), "wb");

        if (fp == nullptr)
        {
            throw std::runtime_error("Cannot open file: " + file_name);
        }

        int rows = static_cast<int>(data.size());
        fwrite(&rows, sizeof(int), 1, fp);

        for (int i = 0; i < rows; ++i)
        {
            int cols = static_cast<int>(data[i].size());
            fwrite(&cols, sizeof(int), 1, fp);

            for (int j = 0; j < cols; ++j)
            {
                int first = data[i][j].first;
                int second = data[i][j].second;
                fwrite(&first, sizeof(int), 1, fp);
                fwrite(&second, sizeof(int), 1, fp);
            }
        }

        fclose(fp);
    }

    //保存top-level-graph节点和对应节点的邻居索引

    void save_top_level_points_to_file(const std::vector<int>& top_level_points, const std::string& file_name)
    {
        FILE *fp=fopen(file_name.c_str(), "wb");

        if(fp==NULL)
        {
            throw std::runtime_error("Cannot open file: "+file_name);
        }

        int num_points=static_cast<int>(top_level_points.size());
        fwrite(&num_points, sizeof(int), 1, fp);
        fwrite(top_level_points.data(), sizeof(int), num_points, fp);

        fclose(fp);
    }

    void save_top_level_indices(const std::vector<std::pair<int, std::vector<int>>>& top_level_indices, const std::string& file_name)
    {
        FILE *fp=fopen(file_name.c_str(),"wb");
        
        if(fp==NULL)
        {
            throw std::runtime_error("Cannot open file: "+file_name);
        }

        int num_pairs=static_cast<int>(top_level_indices.size());
        fwrite(&num_pairs,sizeof(int),1,fp);

        for(const auto& pair:top_level_indices)
        {
            fwrite(&pair.first,sizeof(int),1,fp);
            int num_elements=static_cast<int>(pair.second.size());
            fwrite(&num_elements,sizeof(int),1,fp);
            fwrite(pair.second.data(),sizeof(int),num_elements,fp);
        }
        
        fclose(fp);
    }

    //保存计算的dist[][]
    void save_dist_to_file(const std::vector<std::vector<double>>& data,const std::string& file_name)
    {
        FILE *fp=fopen(file_name.c_str(),"wb");

        if(fp==NULL)
        {
            throw std::runtime_error("Cannot open file: "+file_name);
        } 

        int rows=static_cast<int>(data.size());
        int cols=(rows>0)? static_cast<int>(data[0].size()):0;

        fwrite(&rows,sizeof(int),1,fp);
        fwrite(&cols,sizeof(int),1,fp);
        
        for(int i = 0; i < rows; ++i)
        {
            if(static_cast<int>(data[i].size()) != cols)
            {
                fclose(fp);
                throw std::runtime_error("Mismatch in column sizes in row " + std::to_string(i));
            }
            fwrite(data[i].data(), sizeof(double), cols, fp);
        }
        
        fclose(fp);
    }


    template <typename T>
    std::vector<std::vector<T>> load_from_file(const std::string& file_name)
    {
        std::cout << "load_from_file " << file_name << std::endl;

        FILE *fp = fopen(file_name.c_str(), "rb");

        if (fp == NULL)
        {
            throw GnnsException(std::string("File ") + file_name + "is NOT exist!\n");
        }

        //read the vector dimension
              
        std::vector<std::vector<T>> matrix;
        

        //allocate memory for the vectors
        try {
            while(true)
            {
                int dim;
                size_t read_count=fread(&dim,4,1,fp);
                if(read_count!=1)
                {
                    if(feof(fp)) break;
                    fclose(fp);
                    throw GnnsException("Error reading file "+file_name+"\n");
                }

                std::vector<T> row(dim);

                read_count=fread(row.data(),sizeof(T),dim,fp);
                if(read_count!=static_cast<size_t>(dim))
                {
                    fclose(fp);
                    throw GnnsException("Error reading file "+file_name+"\n");
                }
                matrix.push_back(std::move(row));
            }
            fclose(fp);
        }catch(const std::bad_alloc &ba)
        {
            fclose(fp);
            throw ba;
        }
        return matrix;
    }

    std::vector<std::vector<std::pair<int,int>>> load_start_end_positions_from_file(const std::string& file_name)
    {
        std::vector<std::vector<std::pair<int,int>>> data;

        FILE *fp = fopen(file_name.c_str(), "rb");
        if (fp == NULL)
        {
            throw std::runtime_error("File " + file_name + " does not exist or cannot be opened.");
        }

        try
        {
            int rows, cols;
            size_t read_count = fread(&rows, sizeof(int), 1, fp);
            if (read_count != 1)
            {
                fclose(fp);
                throw std::runtime_error("Error reading file " + file_name);
            }

            data.resize(rows);

            for (int i = 0; i < rows; ++i)
            {
                read_count = fread(&cols, sizeof(int), 1, fp);
                if (read_count != 1)
                {
                    fclose(fp);
                    throw std::runtime_error("Error reading file " + file_name);
                }

                data[i].resize(cols);

                for (int j = 0; j < cols; ++j)
                {
                    int first, second;
                    read_count = fread(&first, sizeof(int), 1, fp);
                    if (read_count != 1)
                    {
                        fclose(fp);
                        throw std::runtime_error("Error reading file " + file_name);
                    }
                    read_count = fread(&second, sizeof(int), 1, fp);
                    if (read_count != 1)
                    {
                        fclose(fp);
                        throw std::runtime_error("Error reading file " + file_name);
                    }

                    data[i][j] = std::make_pair(first, second);
                }
            }

            fclose(fp);
        }
        catch (const std::exception& e)
        {
            fclose(fp);
            throw e;
        }

        return data;
    }

    std::vector<std::vector<std::pair<int,int>>> load_top_level_start_end_positions_from_file(const std::string& file_name)
    {
        std::vector<std::vector<std::pair<int,int>>> data;

        FILE *fp = fopen(file_name.c_str(), "rb");
        if (fp == NULL)
        {
            throw std::runtime_error("File " + file_name + " does not exist or cannot be opened.");
        }

        try
        {
            int rows, cols;
            size_t read_count = fread(&rows, sizeof(int), 1, fp);
            if (read_count != 1)
            {
                fclose(fp);
                throw std::runtime_error("Error reading file " + file_name);
            }

            data.resize(rows);

            for (int i = 0; i < rows; ++i)
            {
                read_count = fread(&cols, sizeof(int), 1, fp);
                if (read_count != 1)
                {
                    fclose(fp);
                    throw std::runtime_error("Error reading file " + file_name);
                }

                data[i].resize(cols);

                for (int j = 0; j < cols; ++j)
                {
                    int first, second;
                    read_count = fread(&first, sizeof(int), 1, fp);
                    if (read_count != 1)
                    {
                        fclose(fp);
                        throw std::runtime_error("Error reading file " + file_name);
                    }
                    read_count = fread(&second, sizeof(int), 1, fp);
                    if (read_count != 1)
                    {
                        fclose(fp);
                        throw std::runtime_error("Error reading file " + file_name);
                    }

                    data[i][j] = std::make_pair(first, second);
                }
            }

            fclose(fp);
        }
        catch (const std::exception& e)
        {
            fclose(fp);
            throw e;
        }

        return data;
    }
    //下载top-level-graph节点和对应节点的邻居索引

    std::vector<int> load_top_level_points_from_file(const std::string& file_name)
    {
        std::vector<int> top_level_points;

        FILE *fp=fopen(file_name.c_str(),"rb");
        if(fp==NULL)
        {
            throw std::runtime_error("File "+file_name+" does not exist or cannot be opened.");
        }

        try
        {
            int num_points;
            size_t read_count = fread(&num_points, sizeof(int), 1, fp);
            if (read_count != 1)
            {
                fclose(fp);
                throw std::runtime_error("Error reading file " + file_name);
            }

            top_level_points.resize(num_points);
            read_count = fread(top_level_points.data(), sizeof(int), num_points, fp);
            if (read_count != static_cast<size_t>(num_points))
            {
                fclose(fp);
                throw std::runtime_error("Error reading file " + file_name);
            }

            fclose(fp);
        }
        catch (const std::exception& e)
        {
            fclose(fp);
            throw e;
        }

        return top_level_points;
    }

    std::vector<std::pair<int,std::vector<int>>> load_top_level_indices_from_file(const std::string& file_name)
    {
        std::vector<std::pair<int, std::vector<int>>> top_level_indices;

        FILE *fp = fopen(file_name.c_str(), "rb");
        if (fp == NULL)
        {
            throw std::runtime_error("File " + file_name + " does not exist or cannot be opened.");
        }

        try
        {
            int num_pairs;
            size_t read_count = fread(&num_pairs, sizeof(int), 1, fp);
            if (read_count != 1)
            {
                fclose(fp);
                throw std::runtime_error("Error reading file " + file_name);
            }

            for (int i = 0; i < num_pairs; ++i)
            {
                int first;
                read_count = fread(&first, sizeof(int), 1, fp);
                if (read_count != 1)
                {
                    fclose(fp);
                    throw std::runtime_error("Error reading file " + file_name);
                }

                int num_elements;
                read_count = fread(&num_elements, sizeof(int), 1, fp);
                if (read_count != 1)
                {
                    fclose(fp);
                    throw std::runtime_error("Error reading file " + file_name);
                }

                std::vector<int> second(num_elements);
                read_count = fread(second.data(), sizeof(int), num_elements, fp);
                if (read_count != static_cast<size_t>(num_elements))
                {
                    fclose(fp);
                    throw std::runtime_error("Error reading file " + file_name);
                }

                top_level_indices.emplace_back(first, std::move(second));
            }

            fclose(fp);
        }
        catch (const std::exception& e)
        {
            fclose(fp);
            throw e;
        }

        return top_level_indices;
    }

    //保存计算的dist[][]
    std::vector<std::vector<double>> load_2d_double_vector_from_file(const std::string& file_name)
    {
        std::vector<std::vector<double>> data;

        FILE *fp = fopen(file_name.c_str(), "rb");
        if (fp == NULL)
        {
            throw std::runtime_error("File " + file_name + " does not exist or cannot be opened.");
        }

        try
        {
            int rows, cols;
            size_t read_count = fread(&rows, sizeof(int), 1, fp);
            if (read_count != 1)
            {
                fclose(fp);
                throw std::runtime_error("Error reading file " + file_name);
            }

            read_count = fread(&cols, sizeof(int), 1, fp);
            if (read_count != 1)
            {
                fclose(fp);
                throw std::runtime_error("Error reading file " + file_name);
            }

            data.resize(rows, std::vector<double>(cols));

            for (int i = 0; i < rows; ++i)
            {
                read_count = fread(data[i].data(), sizeof(double), cols, fp);
                if (read_count != static_cast<size_t>(cols))
                {
                    fclose(fp);
                    throw std::runtime_error("Error reading file " + file_name);
                }
            }

            fclose(fp);
        }
        catch (const std::exception& e)
        {
            fclose(fp);
            throw e;
        }

        return data;
    }

//}

#endif //GNNS_IO_H
