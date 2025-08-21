pipeline {
  agent {
    docker { image 'node:22.18.0-alpine3.22' }
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