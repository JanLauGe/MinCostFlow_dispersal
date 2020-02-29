context("test-funs_constraint_matrix")

data('phillips_problem')

phillips_constraints_flow2 <- create_constraints_matrix(
  edges = phillips_problem,
  total_flow = 2)

phillips_constraints_flow2_lhs <- structure(
  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1),
  .Dim = c(32L, 14L),
  .Dimnames = list(
    c("1", "2",
      "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
      "2", "3", "4", "5", "6", "7", "8", "9", "2", "3", "4", "5", "6",
      "7", "8", "9", "source", "target"),
    c("1", "2", "3", "4", "5",
      "6", "7", "8", "9", "10", "11", "12", "13", "14")))

phillips_constraints_flow2_dir <- c(
  "<=", "<=", "<=", "<=", "<=", "<=", "<=", "<=", "<=", "<=",
  "<=", "<=", "<=", "<=", "==", "==", "==", "==", "==", "==", "==",
  "==", "<=", "<=", "<=", "<=", "<=", "<=", "<=", "<=", "==", "==")

phillips_constraints_flow2_rhs <- c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 1, 1, 1, 1, 1, 1, 1, -2, 2)


test_that('phillips problem data is as expected', {
  # class
  expect_is(phillips_problem, class = 'data.frame')
  # shape
  expect_equal(dim(phillips_problem), c(14, 5))
  # content
  expect_identical(
    data.frame(
      edge_id = seq(1, 14, 1),
      node_from = c(1, 1, 1, 1, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9),
      node_to = c(2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 10, 10, 10, 10),
      edge_capacity = rep(1, 14),
      edge_cost = rep(1, 14)),
    phillips_problem
  )
})

test_that('create_constraints_matrix runs without errors', {
  expect_silent(
    out <- create_constraints_matrix(
      edges = phillips_problem,
      total_flow = 2))
})

test_that('create_constraints_matrix format is correct', {
  expect_equal(
    create_constraints_matrix(
      edges = phillips_problem,
      total_flow = 2) %>%
      length(), 3)
  expect_equal(
    create_constraints_matrix(
      edges = phillips_problem,
      total_flow = 2) %>%
      names(), c('lhs', 'dir', 'rhs'))
})

test_that('create_constraints_matrix for phillips problem yields expected result', {
  expect_equal(
    create_constraints_matrix(edges = phillips_problem, total_flow = 2)[['lhs']],
    phillips_constraints_flow2_lhs)
  expect_equal(
    create_constraints_matrix(edges = phillips_problem, total_flow = 2)[['dir']],
    phillips_constraints_flow2_dir)
  expect_equal(
    create_constraints_matrix(edges = phillips_problem, total_flow = 2)[['rhs']],
    phillips_constraints_flow2_rhs)
})

test_that('add_constraints function runs without errors', {
  expect_silent(
    out <- add_constraints(
      constraints = phillips_constraints_flow2,
      new_lhs = matrix(nrow = 3, ncol = 14, data = 0),
      new_dir = rep('==', times = 3),
      new_rhs = rep(1, times = 3))
  )
})

test_that('add_constraints yields expected result', {
  expect_equal(
    add_constraints(
      constraints = phillips_constraints_flow2,
      new_lhs = matrix(nrow = 3, ncol = 14, data = 0),
      new_dir = rep('==', times = 3),
      new_rhs = rep(1, times = 3)),
    list(
      'lhs' = rbind(
        phillips_constraints_flow2[['lhs']],
        matrix(nrow = 3, ncol = 14, data = 0)),
      'dir' = c(
        phillips_constraints_flow2[['dir']],
        rep('==', times = 3)),
      'rhs' = c(
        phillips_constraints_flow2[['rhs']],
        rep(1, times = 3))
    )
  )
})

