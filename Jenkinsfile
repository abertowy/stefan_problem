pipeline {
  agent {
    docker { image 'alpine:latest' }
  }
  stages {
    stage('whoami') {
      steps {
        sh 'whoami'
        sh 'pwd'
      }
    }
  }
}