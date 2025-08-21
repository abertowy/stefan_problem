pipeline {
  agent any
  stages {
    stage('checkout') {
      steps {
        git(url: 'https://github.com/abertowy/stefan_problem.git', branch: 'bazel_build')
      }
    }

    stage('build') {
      steps {
        sh 'bazel build //stefanproblem:stefanproblem_bazel'
      }
    }

  }
}