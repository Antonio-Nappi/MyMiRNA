import { HttpClient, HttpHeaders, HttpParams } from '@angular/common/http';
import { Injectable } from '@angular/core';

@Injectable({
    providedIn: 'root'
})

export class RequestService {

    constructor(private http: HttpClient) {}

    trimmingStep(body: string) {
        const path = 'http://localhost:8080/trimming';
        return this.http.post(
            path,
            body,
            {
                headers: new HttpHeaders({
                    'Content-Type': 'application/json',
                })
            });
    }

    shortStack(body: string) {
        const path = 'http://localhost:8080/shortstack';
        return this.http.post(
            path,
            body,
            {
                headers: new HttpHeaders({
                    'Content-Type': 'application/json',
                })
            });
    }

    mirnaDetection(body: string) {
        const path = 'http://localhost:8080/mirna';
        return this.http.post(
            path,
            body,
            {
                headers: new HttpHeaders({
                    'Content-Type': 'application/json',
                })
            });
    }
}
