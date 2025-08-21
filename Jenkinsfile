pipeline {
  agent {
    docker { image 'alpine:3.14' }
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