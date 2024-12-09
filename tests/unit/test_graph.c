#include <check.h>
#include <stdlib.h>
#include "../headers/graph.h"
#include "../headers/structs.h"
#include "../mocks/graph_mock.c"

START_TEST(test_initialise_graph)
{
    DatasetInfo dataset;
    dataset.num_vectors = 5;
    dataset.datapoints = (DataPoint *)malloc(dataset.num_vectors * sizeof(DataPoint));
    for (int i = 0; i < dataset.num_vectors; i++) {
        dataset.datapoints[i].category = i;
        dataset.datapoints[i].point_index = i;
        for (int j = 0; j < 100; j++) {
            dataset.datapoints[i].vectors[j] = j * 0.1;
        }
    }

    Graph graph = initialise_graph(&dataset, 10);

    ck_assert_int_eq(graph.num_points, dataset.num_vectors);
    for (int i = 0; i < dataset.num_vectors; i++) {
        ck_assert_int_eq(graph.points[i].index, i);
        ck_assert_int_eq(graph.points[i].category, dataset.datapoints[i].category);
        for (int j = 0; j < 100; j++) {
            ck_assert_float_eq_tol(graph.points[i].coordinates[j], dataset.datapoints[i].vectors[j], 0.001);
        }
        ck_assert_int_eq(graph.points[i].edge_count, 0);
    }
    free(dataset.datapoints);
    free_graph(graph);
}
END_TEST

START_TEST(test_addEdge)
{
    Graph graph;
    graph.num_points = 2;
    graph.points = (Point *)malloc(graph.num_points * sizeof(Point));
    graph.points[0].index = 0;
    graph.points[0].category = 0;
    graph.points[0].coordinates = (float *)malloc(100 * sizeof(float));
    for (int j = 0; j < 100; j++) {
        graph.points[0].coordinates[j] = j * 0.1;
    }
    graph.points[0].edge_count = 0;
    graph.points[0].edges = (int *)malloc(10 * sizeof(int));

    graph.points[1].index = 1;
    graph.points[1].category = 1;
    graph.points[1].coordinates = (float *)malloc(100 * sizeof(float));
    for (int j = 0; j < 100; j++) {
        graph.points[1].coordinates[j] = j * 0.1;
    }
    graph.points[1].edge_count = 0;
    graph.points[1].edges = (int *)malloc(10 * sizeof(int));

    addEdge(&graph.points[0], 1);

    ck_assert_int_eq(graph.points[0].edge_count, 1);
    ck_assert_int_eq(graph.points[0].edges[0], 1);

    addEdge(&graph.points[0], 1);

    ck_assert_int_eq(graph.points[0].edge_count, 1);

    free_graph(graph);
}
END_TEST

START_TEST(test_add_to_dynamic_array)
{
    int *array = NULL;
    int size = 0;

    add_to_dynamic_array(&array, &size, 1);

    ck_assert_int_eq(array[0], 1);
    ck_assert_int_eq(size, 1);

    add_to_dynamic_array(&array, &size, 2);

    ck_assert_int_eq(array[0], 1);
    ck_assert_int_eq(size, 2);
    ck_assert_int_eq(array[1], 2);

    add_to_dynamic_array(&array, &size, 3);

    ck_assert_int_eq(array[0], 1);
    ck_assert_int_eq(size, 3);
    ck_assert_int_eq(array[2], 3);

    free(array);
}
END_TEST

START_TEST(test_filtered_Robust_prune)
{
    Graph* graph = create_graph();

    int *array = (int *)malloc(5 * sizeof(int));
    for (int i = 0; i < 5; i++) {
        array[i] = i;
    }
    int array_size = 5;

    filtered_Robust_prune(graph, 1, array, array_size, 1.3, 5);

    ck_assert_int_eq(graph->points[0].edge_count, 5);
    ck_assert_int_eq(graph->points[1].edge_count, 3);
    ck_assert_int_eq(graph->points[2].edge_count, 5);

    free(array);
}
END_TEST

START_TEST(test_sort_array)
{

    DatasetInfo dataset;
    dataset.num_vectors = 10;
    dataset.datapoints = (DataPoint *)malloc(10 * sizeof(DataPoint));

    for (int i = 0; i < 10; i++) {
        dataset.datapoints[i].category = i;
        dataset.datapoints[i].point_index = i;
        for (int j = 0; j < 100; j++) {
            dataset.datapoints[i].vectors[j] = j * i * 0.1;
        }
    }

    Graph graph = create_random_graph(dataset, 100, 10);

    int *array = (int *)malloc(8 * sizeof(int));
    for (int i = 0; i < 8; i++) {
        array[i] = i;
    }
    int array_size = 8;

    float Xq[100];
    for (int i = 0; i < 100; i++) {
        Xq[i] = i * 0.1;
    }

    sort_array(&graph, array, array_size, Xq);

    ck_assert_int_eq(array[3], 3);
    ck_assert_int_eq(array[4], 4);
    ck_assert_int_eq(array[5], 5);
    ck_assert_int_eq(array[6], 6);
    ck_assert_int_eq(array[7], 7);

    free(array);
}
END_TEST

START_TEST(test_sort_array_based_on_dataset)
{
    DatasetInfo dataset;
    dataset.num_vectors = 10;
    dataset.datapoints = (DataPoint *)malloc(10 * sizeof(DataPoint));

    for (int i = 0; i < 10; i++) {
        dataset.datapoints[i].category = i;
        dataset.datapoints[i].point_index = i;
        for (int j = 0; j < 100; j++) {
            dataset.datapoints[i].vectors[j] = j * 0.1;
        }
    }

    int *array = (int *)malloc(5 * sizeof(int));
    for (int i = 0; i < 5; i++) {
        array[i] = i;
    }
    int array_size = 5;

    float Xq[100];
    for (int i = 0; i < 100; i++) {
        Xq[i] = i * 0.1;
    }

    sort_array_based_on_dataset(&dataset, array, array_size, Xq);

    ck_assert_int_eq(array[0], 0);
    ck_assert_int_eq(array[1], 1);
    ck_assert_int_eq(array[2], 2);
    ck_assert_int_eq(array[3], 3);
    ck_assert_int_eq(array[4], 4);

    free(dataset.datapoints);
    free(array);
}
END_TEST

START_TEST(test_edgeExists_empty_graph)
{
    Graph graph;
    graph.num_points = 0;
    graph.points = (Point *)malloc(sizeof(Point));
    graph.points[0].edge_count = 0;
    graph.points[0].index = 0;
    graph.points[0].edges = (int *)malloc(sizeof(int));

    ck_assert_int_eq(edgeExists(&graph.points[0], 1), 0);
}
END_TEST

START_TEST(test_filtered_Robust_prune_empty_array)
{
    Graph *graph = create_graph();

    int *array = NULL;
    int size = 0;

    
    int check = filtered_Robust_prune(graph, 1, array, size, 1.5, 10);

    ck_assert_int_eq(check, 0);


    free_graph(*graph);
}
END_TEST

START_TEST(test_sort_array_empty_array)
{
    Graph graph;
    graph.num_points = 0;
    graph.points = NULL;

    int *array = NULL;
    int size = 0;
    float Xq[100];
    for (int i = 0; i < 100; i++) {
        Xq[i] = i * 0.1;
    }

    sort_array(&graph, array, size, Xq);

    ck_assert_ptr_eq(array, NULL);
    ck_assert_int_eq(size, 0);
}
END_TEST

START_TEST(test_sort_array_based_on_dataset_empty_dataset)
{
    DatasetInfo dataset;
    dataset.num_vectors = 0;
    dataset.datapoints = NULL;

    int *array = (int *)malloc(5 * sizeof(int));
    for (int i = 0; i < 5; i++) {
        array[i] = i;
    }
    float number = 4.0;

    int check = sort_array_based_on_dataset(&dataset, array, 5, &number);

    ck_assert_ptr_eq(dataset.datapoints, NULL);
    ck_assert_int_eq(check, 0);

    free(array);
}
END_TEST

START_TEST(test_get_the_difference)
{
    int a[] = {1, 2, 3};
    int b[] = {2, 3, 6};
    int a_size = 3;
    int b_size = 3;
    int c_size = 0;
    int expected[] = {3};

    int *result = get_the_difference(a, a_size, b, b_size, &c_size);
    printArray(result, c_size);
    ck_assert_int_eq(result[0], 6);
    ck_assert_int_eq(c_size, 4);
    free(result);
}
END_TEST


START_TEST(test_arrayContains)
{
    int array[] = {1, 2, 3};
    int size = 3;
    int value = 2;

    ck_assert_int_eq(arrayContains(array, size, value), 1);
    ck_assert_int_eq(!arrayContains(array, size, 4), 1);
}
END_TEST

START_TEST(test_squared_euclidean_distance)
{
    float a[] = {1.0, 2.0, 3.0};
    float b[] = {4.0, 5.0, 6.0};
    float expected = 27.0;
    double result = squared_euclidean_distance(a, b, 3);

    ck_assert_float_eq(squared_euclidean_distance(a, b, 3), expected);
}
END_TEST

Suite* graph_suite(void)
{
    Suite* s = suite_create("graph_suite");
    TCase* c_core = tcase_create("Core");

    tcase_add_test(c_core, test_initialise_graph);
    tcase_add_test(c_core, test_addEdge);
    tcase_add_test(c_core, test_add_to_dynamic_array);
    tcase_add_test(c_core, test_filtered_Robust_prune);
    tcase_add_test(c_core, test_sort_array);
    tcase_add_test(c_core, test_sort_array_based_on_dataset);
    tcase_add_test(c_core, test_edgeExists_empty_graph);
    tcase_add_test(c_core, test_filtered_Robust_prune_empty_array);
    tcase_add_test(c_core, test_sort_array_empty_array);
    tcase_add_test(c_core, test_sort_array_based_on_dataset_empty_dataset);
    tcase_add_test(c_core, test_get_the_difference);
    tcase_add_test(c_core, test_arrayContains);
    tcase_add_test(c_core, test_squared_euclidean_distance);



    suite_add_tcase(s, c_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite* s = graph_suite();
    SRunner* sr = srunner_create(s);
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
