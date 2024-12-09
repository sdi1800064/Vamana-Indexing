#include <check.h>
#include <stdlib.h>
#include "tests/unit/test_graph.c"

// Include suite functions from individual test files
Suite *graph_suite(void);
// Add more as needed

int main(void) {
    int num_failed;
    SRunner *sr = srunner_create(graph_suite());

    // // Add other test suites
    // srunner_add_suite(sr, file2_suite());

    srunner_run_all(sr, CK_VERBOSE);
    num_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
